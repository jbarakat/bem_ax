import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import fnmatch as fn

# list global trajectory files (number corresponds to number of cells placed)
fsteps = sorted(glob.glob('../output/nodexr*'))
nsteps = len(fsteps)

# import data
for f in fsteps :
	data = np.loadtxt(f,unpack=True)
	x = []
	y = []
	x.append(data[0,:]-data[0,0])
	y.append(data[1,:])
	x = np.column_stack(x)
	y = np.column_stack(y)
	plt.plot(x, y, '-o')

plt.xlabel('x')
plt.ylabel('r')
#plt.xlim((0,0.5))
#plt.ylim((0.4, 0.6))
plt.legend(loc='upper right')
#plt.savefig('plot.png')
plt.show()
