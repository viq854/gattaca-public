import logging
import numpy as np

from scipy.misc import logsumexp
from scipy.special import gammaln, psi
from sklearn import cluster

from misc import AttributeDict, gen_batches

np.seterr(all='warn')
tol=1e-4

default_init = {
  'm0'    : None,
  'm'     : None,
  'W0'    : None,
  'Winv0' : None,
  'W'     : None,
}

# ----------------------------------------------------------------------------
# clustering code

def variational_em(X_all, K, init=default_init, 
                   max_epoch=100, eps=1e-4, n_minibatch=1000, rho=1):
  d, n_all = X_all.shape[0], X_all.shape[1]
  n = n_minibatch
  it, llik, diff = 0, float('-inf'), float('inf')
  # hyperparameters
  hp = AttributeDict()
  hp.beta0  = 1e-3
  hp.nu0    = d+2

  variances = np.var(X_all, axis=1)
  hp.m0     = np.mean(X_all,axis=1) if init['m0'] is None \
              else init['m0']
  hp.W0    = np.diag(1./(variances*d)) if init['W0'] is None else init['W0']
  hp.Winv0 = np.diag(variances*d) if init['Winv0'] is None else init['Winv0']

  logging.info('Initializing clustering model via k-means...')

  # parameters of approximate model q
  pa = AttributeDict()

  # q(mu_k, Lam | nu,W) = q(mu_k | m_k, beta_k * Lam_k) q(Lam | nu,W) (NIW)
  pa.beta  = 1e-3*np.ones(K,)
  pa.nu    = d*np.ones(K,)+2
  pa.W     = 0.1*np.tile(np.eye(d).reshape((1,d,d)), (K,1,1))
  pa.m, pa.W, logpi = _init_means_and_cov(X_all, K, hp)

  if init['m'] is not None: pa.m = init['m']
  if init['W'] is not None: pa.W = init['W']

  # gradients (for step size calculations)
  eg = AttributeDict() # expected gradients
  eh = AttributeDict() # expected grad dot product
  st = AttributeDict() # step sizes
  wi = AttributeDict() # window sizes (tau in Ranganath et al. (2013))
  dv = AttributeDict() # stochastic gradients
  for param in pa.keys():
    eg[param] = None
    st[param] = None
    wi[param] = None
    dv[param] = None
  first_it = True

  asgn = np.empty(n_all,dtype=np.int)

  for epoch in xrange(max_epoch):
    for i, (n, batch) in enumerate(gen_batches(n_all, n_minibatch)):
      # load minibatch
      # X = X_all
      X = X_all[:, batch]
      batch_adj = float(n_all) / n
      logging.info('Epoch: %d/%d. Log-likelihood: %f' % (epoch, max_epoch, llik))

      logp  = np.zeros((n,1))

      # compute q(s | p)

      # compute the expectations

      # # Exp [ ln(Lam) ] and Exp [ ln(pi) ]
      _, logDetW = np.linalg.slogdet(pa.W) # (K,)
      expLogDetLam = _expLogDetLam(pa.nu, logDetW, d)

      # compute sufficient statistics
      # dx = pa.m.reshape((d,K,1)) - X.reshape((d,1,n)) # (d,K,n)
      # dxWdx = np.einsum('ikn, kij, jkn -> kn', dx, pa.W, dx) # (K,n)
      dxWdx = np.empty((K,n), dtype=np.float32)
      for k in xrange(K):
        dx = pa.m[:,k].reshape((d,1)) - X # (d,n)
        dxWdx[k] = np.einsum('in, ij, jn -> n', dx, pa.W[k], dx)
      nudxWdx = pa.nu.reshape((K,1)) * dxWdx

      logp = 0.5*( expLogDetLam.reshape((K,1))
                   + 2*logpi.reshape((K,1))
                   - nudxWdx 
                   - d/pa.beta.reshape((K,1))
                  ) # (K,n)
      logp -= logsumexp(logp, axis=0).reshape(1,n) # (k,n)
      logpk = logsumexp(logp, axis=1) # (k,)
      logp_avg = logp - logpk.reshape(K,1) #(k,n)
      p = np.exp(logp)
      pk = np.exp( logpk ) # (k,)
      p_avg = np.exp(logp_avg)

      # compute q(mu, Lam, pi)

      # helpers
      Nk = batch_adj * pk # (K,)
      # Nk = pk # (K,)
      Xk, Sk = _gaussian_estimates(X, p)

      dxm = Xk - hp.m0.reshape(d,1) #(d,K)
      dxmdxmT = np.einsum('dk,ck->kdc', dxm, dxm)
      beta_wts = (hp.beta0 * Nk) / (hp.beta0 + Nk)

      # gradients
      dv.beta = hp.beta0 + Nk - pa.beta
      dv.m = ( hp.beta0*hp.m0.reshape((d,1)) + Nk*Xk ) / (hp.beta0 + Nk).reshape(1,K) - pa.m
      Winv = hp.Winv0.reshape(1,d,d) + Nk.reshape(K,1,1) * Sk + beta_wts.reshape(K,1,1) * dxmdxmT
      dv.W = np.linalg.inv(Winv) - pa.W
      dv.nu = hp.nu0 + Nk - pa.nu

      # update expected gradients
      if first_it:
        for param in pa.keys():
          eg[param] = dv[param]
          eh[param] = 1.1*np.sum(dv[param]**2) # initial step size of 0.1
          wi[param] = n
        first_it  = False
      else:
        for param in pa.keys():
          tau_inv = 1./wi[param]
          eg[param] = (1-tau_inv) * eg[param] + tau_inv * dv[param]
          eh[param] = (1-tau_inv) * eh[param] + tau_inv * np.sum(dv[param]**2)

      # compute step size, windows size, and update parameters
      for param in pa.keys():
        if eh[param] > 1e-4:
          st[param] = np.sum(eg[param]**2) / eh[param]
        if n >= n_all:
          st[param] = 1.
        # st[param] = 0.8
        wi[param] = wi[param] * (1-st[param]) + 1
        pa[param] = pa[param] + st[param] * dv[param]

      # update assignments
      asgn[batch] = np.argmax(logp,axis=0)

      # llik_new = _full_bound(X_all, pa, hp, K)
      llik_new = _bound(logp, pa, hp, logpi, Nk, Xk, Sk)
      diff = llik_new - llik

      ## m step:

      # compute pi

      logpi = logpk - np.log(n) # (K,)
          
      llik = llik_new
      it += 1


    if abs(diff) < eps:
      break

  _, logDetW = np.linalg.slogdet(pa.W) # (K,)
  expLogDetLam = _expLogDetLam(pa.nu, logDetW, d)

  # compute responsibilities for entire dataset
  for n, batch in gen_batches(n_all, n):
    X = X_all[:,batch]
    dxWdx = np.empty((K,n))
    dxk = np.empty((K,n))
    for k in xrange(K):
      dx = pa.m[:,k].reshape((d,1)) - X # (d,n)
      dxWdx[k] = np.einsum('in, ij, jn -> n', dx, pa.W[k], dx)
    nudxWdx = pa.nu.reshape((K,1)) * dxWdx

    logp = 0.5*( expLogDetLam.reshape((K,1)) 
                 + 2*logpi.reshape((K,1))
                 - nudxWdx - d/pa.beta.reshape((K,1))
                ) # (K,n)
    logp -= logsumexp(logp, axis=0).reshape(1,n)
    asgn[batch] = np.argmax(logp,axis=0)

  Xk, Sk = _gaussian_estimates(X, np.exp(logp))

  Sigma = np.linalg.inv(pa.W)
  Sigma /= (pa.nu.reshape(K,1,1))

  return pa.m, Sigma, asgn, llik

# ----------------------------------------------------------------------------
# initialization

def _init_means_and_cov(X, K, hp):
  d = X.shape[0]
  kmeans = cluster.KMeans(n_clusters=K, verbose=False)
  kmeans.fit(X.T)
  m = kmeans.cluster_centers_.T
  z = kmeans.labels_
  W = np.empty((K,d,d))
  zeros = np.zeros((d,d))
  pi = np.empty(K,)
  for k in range(K):
    Xk = X[:,z==k]
    nk = np.count_nonzero(z==k)
    if nk > 1:
      pi[k] = float(nk)
      cov = np.cov(Xk)
      mk = np.mean(X[:,z==k], axis=1)
      dmk = mk - hp.m0
      beta_wt = (hp.beta0 * nk) / (hp.beta0 + nk)
      Winv  = hp.Winv0 + nk * cov + beta_wt * np.einsum('i,j->ij', dmk, dmk)
      W[k] = np.linalg.inv(Winv)
    else:
      pi[k] = 1e-4
      W[k] = hp.W0
  pi /= np.sum(pi)
  return m, W, np.log(pi)

def _init_means(X, K, eps=7.5):
  m = cluster.KMeans(n_clusters=K).fit(X.T).cluster_centers_.T
  m += eps*np.random.randn(*m.shape)
  return m

def _init_grad(X):
  gr = AttributeDict()

  eg.m     = (X_all.mean(axis=1).reshape(d,1) + 10*np.random.randn(d,K)) \
              if m_init is None else m_init
  pa.m     = _init_means(X_all,K) if m_init is None else m_init
  pa.beta  = np.random.rand(K,)
  pa.nu    = 2*d*np.ones(K,)+2
  pa.W     = np.tile(np.eye(d).reshape((1,d,d)), (K,1,1))
  # W     = _init_W(X,K)

# ----------------------------------------------------------------------------
# helpers

def _bound(logp, pa, hp, logpi, Nk, Xk, Sk):
  d, K = Xk.shape
  pk = Nk
  
  _, logDetW = np.linalg.slogdet(pa.W)
  expLogDetLam = _expLogDetLam(pa.nu, logDetW, d)

  # Exp[ log p(X | Z,mu,Lam) ]
  dxmk = Xk - pa.m # (d,K)
  dxmkWdxmk = np.einsum('ik, kij, jk -> k', dxmk, pa.W, dxmk) # (K,)
  tmp = expLogDetLam - d / pa.beta \
      - pa.nu * np.trace(
          np.einsum('kab,kbc->kac', Sk, pa.W), 
          axis1=1, axis2=2
        ) \
      - pa.nu * dxmkWdxmk \
      - d * np.log(2*3.14)
  expLogpX = 0.5 * np.sum(Nk * tmp)

  # Exp[ log p(mu,Lam) ]
  dm = pa.m - hp.m0.reshape(d,1) # (d,K)
  dmWdm = np.einsum('ik, kij, jk -> k', dm, pa.W, dm) # (K,)

  # Wishart log-normalizing constant
  _, logDetW0 = np.linalg.slogdet(hp.W0)
  logB = _log_wishart_Z(hp.nu0, logDetW0, d)

  trWinv0Wk = np.trace(
                np.einsum('ab,kbc->kac', hp.Winv0, pa.W), 
                axis1=1, axis2=2
              ) # (K,)
  
  tmp = 0.5 * ( d * np.log(hp.beta0 / (2*3.14))
                + expLogDetLam
                - d*hp.beta0 / pa.beta
                - hp.beta0 * pa.nu * dmWdm
              ) \
      + logB \
      + 0.5*(hp.nu0 - d - 1) * expLogDetLam \
      - 0.5 * pa.nu * trWinv0Wk
  expLogpMuLam = np.sum(tmp, axis=0)

  # Exp[ log p(Z | pi) ]
  expLogpZ = np.sum(pk * logpi, axis=0)

  # Exp[ log q(Z) ]
  expLogqZ = np.sum(np.exp(logp) * logp)

  # Exp[ log q(mu, Lam) ]

  # log-normalizing constants B(W_k, nu_k)
  logB = _log_wishart_Z(pa.nu, logDetW, d, K)

  # entropy of Wishart; dim = (K,)
  entLam = - logB \
           - 0.5 * (pa.nu - d - 1) * expLogDetLam \
           + 0.5 * pa.nu * d

  tmp = 0.5 * ( expLogDetLam + d * np.log(pa.beta/(2*3.14)) - d ) \
      - entLam
  expLogqMuLam = np.sum(tmp, axis=0)

  # compute bound
  llik_new = expLogpX + expLogpZ + expLogpMuLam \
           - expLogqZ - expLogqMuLam
  return llik_new

# ----------------------------------------------------------------------------

def _gaussian_estimates(X, p):
  d, n = X.shape
  K, n = p.shape

  Xk = np.empty((d,K))
  Sk = np.empty((K,d,d))
  for k in xrange(K):
    sum_p = np.sum(p[k])
    if sum_p:
      p_avg = p[k] / sum_p
      Xpk = X * p_avg # (d,n)
      Xk[:,k] = np.sum(Xpk,axis=1)

      dx = X - Xk[:,k].reshape(d,1) # (d,n)
      dxp = dx * p_avg.reshape(1,n)  
      Sk[k] = np.einsum('dn,cn->dc', dxp, dx)
    else:
      Xk[:,k] = np.zeros(d,)
      sigma = 1e-4 * np.eye(d)

  return Xk, Sk

def _expLogDetLam(nu, logDetW, d):
  """Expectation of log(det(Lam)) under q(mu, Lam | m, beta, W)."""
  return np.array([
             np.sum([ psi( 0.5*(n+1-(di+1)) ) for di in range(d) ])
             for n in nu
           ]) \
         + d * np.log(2) \
         + logDetW # (K,)

def _log_wishart_Z(nu, logDetW, d, K=None):
  """Log normalizing constant for Wishart distribution."""
  if not K:
    return - nu/2. * logDetW \
           - nu*d/2. * np.log(2.) \
           - d*(d-1)/4. * np.log(3.14) \
           - np.sum([ gammaln( 0.5*(nu+1-(di+1)) ) for di in range(int(d)) ])
  else:
    return - nu/2. * logDetW \
           - nu*d/2. * np.log(2.) \
           - d*(d-1)/4. * np.log(3.14) \
           - np.array([
              np.sum([ gammaln( 0.5*(nu[k]+1-(di+1)) ) for di in range(int(d)) ])
              for k in range(K)
             ]) # (K,)