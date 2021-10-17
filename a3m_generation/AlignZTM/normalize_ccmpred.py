import sys, math
import numpy as np

f = open(sys.argv[1],'r')
d = np.array(f.read().strip().split(),np.float)
f.close()
# with open() as f: d = np.array(f.read().strip().split(),np.float)

dim = int(math.sqrt(len(d)))
assert dim**2 == len(d)
d -= d.mean()
d /= d.std()
d = d.reshape((dim,dim))

f = open(sys.argv[2],'w')
f.write('\n'.join(['\t'.join([str(d[i,j]) for j in range(dim)]) for i in range(dim)]) + '\n')
f.close()
# with open(sys.argv[2],'w') as f: f.write('\n'.join(['\t'.join([str(d[i,j]) for j in range(dim)]) for i in range(dim)]) + '\n')