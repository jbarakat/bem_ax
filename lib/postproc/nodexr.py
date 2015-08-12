import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import fnmatch as fn

# declare arrays
x = []
y = []

# list global trajectory files (number corresponds to number of cells placed)
fsteps = sorted(glob.glob('../output/nodexr*'))
nsteps = len(fsteps)

# import data
for f in fsteps :
  data = np.loadtxt(f,unpack=True)
  x.append(data[0,:])
  y.append(data[1,:])

x = np.column_stack(x)
y = np.column_stack(y)

# plot profiles
for i in range(nsteps) :
	#plt.plot(x[:,i],y[:,i],'-k',label='ts = '+str(i))
	plt.plot(x[:,i],y[:,i],'-o',label='ts = '+str(i))


plt.xlabel('x')
plt.ylabel('r')
#plt.xlim((0,150))
#plt.ylim((0,15))
plt.legend(loc='upper right')
#plt.savefig('plot_dx.png')
plt.show()

## compute inter-particle distance
#d = np.array([ [ x[i,j] - x[i,j-1] for i in range(len(x)) ] for j in range(nsteps) ])
#d = d.T
#
#for i in range(len(x)) :
#  for j in [0] :
#    d[i,j] = d[i,j] + 15
#
#d = d*(float(nsteps)/15.0)
#
## compute cell position  relative to cell 1
#dx = np.array([ [ x[i,j] - x[i,0] for i in range(len(x)) ] for j in range(nsteps) ])
#dx = dx.T
#
## compute the initial position
#dx0 = np.array([ x[0,j] - x[0,0] for j in range(nsteps)])
#dx0 = dx0.T
#dx0 = np.tile(dx0,(2,1))
#t0 = np.array([0,150])
#t0 = t0.T
##print dx0
##print t0
#
## plot distance relative to cell 1
#for i in range(nsteps) :
#  plt.plot(t[:,i],dx[:,i],label='cell '+str(i+1))
#  plt.plot(t0,dx0[:,i],':k')
#
#plt.xlabel('time')
#plt.ylabel('COM position (relative to cell 1)')
#plt.xlim((0,150))
#plt.ylim((0,15))
#plt.legend(loc='upper right')
#plt.savefig('plot_dx.png')
##plt.show()
#
### plot inter-particle distance
##for i in range(nsteps) :
##  plt.plot(t[:,i],d[:,i],label=str(i+1)+'-'+str(i+2))
##
##plt.xlabel('time, Ut/(2R+<d>)')
##plt.ylabel('center-to-center distance, d/R')
##plt.legend(loc='upper right')
##plt.savefig('plot_dx.png')
###plt.show()
#
#
## plot global trajectories
## plt.plot(t,x)
## plt.xlabel('time')
## plt.ylabel('position')
## plt.savefig('plot_x.png')
#
