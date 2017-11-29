#### Configurable settings

USE_OPENMP=true
BITS_PER_COUNTER=8
HAVE_ZLIB=1

#####
CMPH_DIR=../../third-party/cmph
CMPH_DIR_ABS=$(realpath $(CMPH_DIR))
BF_DIR=../../third-party/libbf
BF_DIR_ABS=$(realpath $(BF_DIR))
SEQAN_INCLUDE_DIR=../../third-party/seqan/include/

PROG=           gattaca
PROG_MINHASH=   gattaca-minh

all: deps_install.done $(PROG) $(PROG_MINHASH) toplevel

deps_install.done:
	@echo "Installing the CMPH library..."
	@cd $(CMPH_DIR) && \
	chmod u+x configure && \
	./configure --prefix=$(CMPH_DIR_ABS) && \
	make && make install
	@echo "CMPH install done."
	@echo "Installing the bf library..."
	@cd $(BF_DIR) && \
        chmod u+x configure && \
        ./configure --prefix=$(BF_DIR_ABS) && \
        make && make install
	@echo "LIBBF install done."
	@touch $@

CC=		g++
CFLAGS=		-msse4.1 -o rand_sse -Wall -O3 -std=c++0x
DFLAGS=		-DSEQAN_HAS_ZLIB=$(HAVE_ZLIB) -DBITS_PER_COUNTER=$(BITS_PER_COUNTER)

ifeq ($(USE_OPENMP), true)
CFLAGS   +=-fopenmp
endif

DEPS=		seq/types.h seq/city.h seq/io.h	
DEPS_INDEX=	index/index.h index/count_min_sketch.h index/mph.h index/basic_map.h index/minhash_sketch.h
DEPS_SELECTION=	selection/minhash.h selection/minhash_params.h

SOURCES_SEQ=	seq/io.cc \
		seq/city.cc \

SOURCES_INDEX=	index/main.cc \
		index/index.cc \

SOURCES_SELECTION= selection/main.cc
		
OBJDIR=		obj
OBJDIR_SEQ=     $(OBJDIR)/seq
OBJDIR_INDEX= $(OBJDIR)/index
OBJDIR_SELECTION= $(OBJDIR)/selection
_OBJS_SEQ=	$(SOURCES_SEQ:.cc=.o)
_OBJS_INDEX=	$(SOURCES_INDEX:.cc=.o)
_OBJS_SELECTION=    $(SOURCES_SELECTION:.cc=.o)
OBJS_SEQ= 	$(patsubst %,$(OBJDIR)/%,$(_OBJS_SEQ))
OBJS_INDEX=$(patsubst %,$(OBJDIR)/%,$(_OBJS_INDEX))
OBJS_SELECTION=$(patsubst %,$(OBJDIR)/%,$(_OBJS_SELECTION))
OBJDIR_GUARD = 		@mkdir -p $(OBJDIR)
OBJDIR_GUARD_SEQ = 	@mkdir -p $(OBJDIR_SEQ)
OBJDIR_GUARD_INDEX = 	@mkdir -p $(OBJDIR_INDEX)
OBJDIR_GUARD_SELECTION =    @mkdir -p $(OBJDIR_SELECTION)

INCLUDES=	-I/usr/local/include -I$(SEQAN_INCLUDE_DIR) -I$(CMPH_DIR)/include -I$(BF_DIR)/include
LIBS=		-lrt -lm -L$(CMPH_DIR)/lib -lcmph -L$(BF_DIR)/lib -lbf -lz -Wl,-R$(BF_DIR_ABS)/lib -Wl,-R$(CMPH_DIR_ABS)/lib
LOG=0

$(OBJDIR_SEQ)/%.o: seq/%.cc $(DEPS)
		$(OBJDIR_GUARD_SEQ)
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

$(OBJDIR_INDEX)/%.o: index/%.cc $(DEPS) $(DEPS_INDEX)
		$(OBJDIR_GUARD_INDEX)
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

$(OBJDIR_SELECTION)/%.o: selection/%.cc $(DEPS) $(DEPS_SELECTION)
		$(OBJDIR_GUARD_SELECTION)
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

$(PROG): bits.check $(OBJS_SEQ) $(OBJS_INDEX)
		@mkdir -p $(OBJDIR)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS_SEQ) $(OBJS_INDEX) -o $@ $(LIBS)

$(PROG_MINHASH): bits.check $(OBJS_SEQ) $(OBJS_SELECTION)
		@mkdir -p $(OBJDIR)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS_SEQ) $(OBJS_SELECTION) -o $@ $(LIBS)

toplevel: $(PROG) $(PROG_MINHASH)
	cp $(PROG) ../../.
	cp $(PROG_MINHASH) ../../.	

debug: DFLAGS += -g -DDEBUG -DLOG=$(LOG)
debug: $(PROG)

bits.check:
ifeq ($(shell test $(BITS_PER_COUNTER) -gt 32; echo $$?),0)
		BITS_PER_COUNTER=64
		@echo WARNING:Counters requiring more than 32 bits will default to 64 bits.
endif
	
clean:
		rm -f $(OBJDIR)/*.o $(OBJDIR_SEQ)/*.o $(OBJDIR_INDEX)/*.o $(OBJDIR_SELECTION)/*.o $(PROG) $(PROG)-* *~ *.a *.done
