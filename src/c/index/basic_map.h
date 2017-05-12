#ifndef BASICMAP_H_
#define BASICMAP_H_

#include "../seq/types.h"
#include <unordered_map>

// multi-stream exact counts table based on a basic hash map
class basic_table_t : public counts_table_t {
public:
	std::unordered_map<kmer_2bit_t, counter_t*> map;
	int n_streams;
	int n_dropped;

	basic_table_t(): basic_table_t(1) {}

	basic_table_t(const uint32 n_input_streams) {
		n_streams  = n_input_streams;
		n_dropped = 0;
	}
	
	// ---- interface ---- //
	virtual int get_n_streams() { 
		return n_streams;
	}

	virtual void insert(const kmer_2bit_t& key, const int stream_id)  { 
		std::unordered_map<kmer_2bit_t, counter_t*>::iterator entry = map.find(key);
		if(entry == map.end()) { // insert new entry
			counter_t* stream_counts = new counter_t[n_streams];
			stream_counts[stream_id] = 1;
			map[key] = stream_counts;
		} else { // update counter
			if(entry->second[stream_id] + 1ULL < (1ULL << BITS_PER_COUNTER)) {
				entry->second[stream_id]++;
			} else {
				n_dropped++;
			}
		}
	}

	virtual counter_t lookup(const kmer_2bit_t& key, const int stream_id) const {
		 std::unordered_map<kmer_2bit_t, counter_t*>::const_iterator entry = map.find(key);
		if(entry == map.end()) { 
			return 0; 
		} else {
			std::cout << key << " " << entry->second[stream_id] << "\n";
			return entry->second[stream_id];
		}
	}
	
	// return the kmer keys that occur more than 'count' times in at least one sample
	virtual  void get_freq_keys(std::vector<kmer_2bit_t>& keys, const int count) {
		for (std::unordered_map<kmer_2bit_t, counter_t*>::const_iterator it = map.begin(); it != map.end(); ++it) {
			kmer_2bit_t key = it->first;
			for(int i = 0; i < n_streams; i++) {
				if(it->second[i] > count) {
					keys.push_back(key);
					break;
				}
			}
		}
	}

	void get_key_values(std::vector<kmer_2bit_t>& keys, std::vector<counter_t>& key_counts) {
		for (std::unordered_map<kmer_2bit_t, counter_t*>::const_iterator it = map.begin(); it != map.end(); ++it) {
			kmer_2bit_t key = it->first;
			if(it->second[0] >= 1) {
				 keys.push_back(key);
				 key_counts.push_back(it->second[0]);
			}
		}
	}

	// return the kmer keys that occur more than 'count' times without loading the full table
	static void get_freq_keys_from_file_per_sample(const std::string& idx_fname, std::vector<std::vector<kmer_2bit_t>>& keys, const int count) {
		std::ifstream fileIDX;
		fileIDX.open(idx_fname.c_str(), std::ios::in | std::ios::binary);
		if (!fileIDX.is_open()) {
			std::cerr << "ERROR: Could not open file: " << idx_fname << "\n";
			exit(1);
		}
		std::cout << "Loading..." << idx_fname << "\n";
		int n_streams;
		fileIDX.read(reinterpret_cast<char*>(&n_streams), sizeof(n_streams));
		long long int map_size;
		fileIDX.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
		std::cout << "Table size: " << map_size << "\n";
		std::cout << "Streams: " << n_streams << "\n";			

		keys.resize(n_streams);	
		while(map_size > 0) {
			kmer_2bit_t key;
			fileIDX.read(reinterpret_cast<char*>(&key), sizeof(key));
			counter_t c;
			for(int n = 0; n < n_streams; n++) {
				fileIDX.read(reinterpret_cast<char*>(&c), sizeof(c));
				if(c > count) {
					keys[n].push_back(key);
				}
			}
			map_size--;
		}
		fileIDX.close();
	}

	virtual ~basic_table_t() {
		for (std::unordered_map<kmer_2bit_t, counter_t*>::const_iterator it = map.begin(); it != map.end(); ++it) {
			free(it->second);
		}
	}

	virtual void clear() {}

	// write the table to file
	virtual void save_to_file(const std::string& fname, const int n_count_bits) {
		std::ofstream file;
		file.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::app);
		file.write(reinterpret_cast<char*>(&n_streams), sizeof(n_streams));
		long long int map_size = map.size();
		file.write(reinterpret_cast<char*>(&map_size), sizeof(map_size));
		for (std::unordered_map<kmer_2bit_t, counter_t*>::const_iterator it = map.begin(); it != map.end(); ++it) {
			file.write(reinterpret_cast<const char*>(&(it->first)), sizeof(it->first));
			file.write(reinterpret_cast<const char*>(&(it->second[0])), n_streams*sizeof(it->second[0]));
		}
		file.close();
	}
	
	// load the table from file
	virtual long int load_from_file(const std::string& fname, long int file_offset) {
		std::ifstream file;
		file.open(fname.c_str(), std::ios::in | std::ios::binary);
		file.seekg(file_offset, file.beg);
		file.read(reinterpret_cast<char*>(&n_streams), sizeof(n_streams));
		long long int map_size;
		file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
		std::cout << "Map size " << map_size << "\n";
		while(map_size > 0) {
			kmer_2bit_t key;
			file.read(reinterpret_cast<char*>(&key), sizeof(key));
			counter_t* stream_counts = new counter_t[n_streams];
			file.read(reinterpret_cast<char*>(&stream_counts[0]), n_streams*sizeof(stream_counts[0]));
			map[key] = stream_counts;
			map_size--;
		}
		long int s = file.tellg();
		file.close();
		return s;
	}
	
	virtual void print_stats() {
		std::cout << "Table size: " << map.size() << ", dropped: " << n_dropped <<"\n";
	}
};

#endif
