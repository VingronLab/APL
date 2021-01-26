import numpy #as np
#import numpy.random
import torch

def svd_torch(x):
	x2 = torch.from_numpy(numpy.array(x))
	u, s, v = torch.svd(x2)
	u1 = u.numpy()
	s1 = s.numpy()
	v1 = v.numpy()
	return u1,s1,v1
