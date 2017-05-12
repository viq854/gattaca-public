import numpy as np

# ----------------------------------------------------------------------------
# object to hold parameters

class AttributeDict(dict):
  def __getattr__(self, attr):
    if attr in self.__dict__:
      return self.attr
    elif attr in self:
      return self[attr]
    else:
      raise ValueError('Invalid attribute')

  def __setattr__(self, attr, value):
    if attr in self.__dict__:
      self.attr = value
    else:
      self[attr] = value

# ----------------------------------------------------------------------------
# batch manipulations

def gen_batches(n_data, n_batchsize):
  range_v = range(n_data)
  np.random.shuffle(range_v)
  n_batches = int(np.ceil( float(n_data) / n_batchsize ))
  for i in xrange(n_batches):
    batch = range_v[i*(n_batchsize):(i+1)*(n_batchsize)]
    yield len(batch), batch