import numpy #as np
#import numpy.random
import torch
from scipy.sparse.linalg import svds
from numpy import array

def svd_torch(x):
  x = x.todense()
  x2 = torch.from_numpy(numpy.array(x))
  u, s, v = torch.svd(x2)
  u1 = u.numpy()
  s1 = s.numpy()
  v1 = v.numpy()
  return u1,s1,v1


def svd_linalg_torch(x):
  x2 = torch.from_numpy(numpy.array(x))
  u, s, v = torch.linalg.svd(x2, full_matrices=False, out=None)
  u1 = u.numpy()
  s1 = s.numpy()
  v1 = v.numpy()
  return u1,s1,v1

def svds_scipy(x, k=6,which = 'LM',ncv=None, tol=0, v0=None, maxiter=None, return_singular_vectors=True, solver='propack'):
	u, s, vh = svds(x, k = int(k), which = which, ncv=ncv, tol=tol, v0=v0, maxiter=maxiter, return_singular_vectors=return_singular_vectors, solver=solver)
	v = vh.transpose()
	return u, s, v
