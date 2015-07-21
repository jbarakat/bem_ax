import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import fnmatch as fn

# declare variables
fn = './str.dat'
X = []
Y = []

# import data
data = np.loadtxt(fn, unpack = True, skiprows = 1)
X.append(data[0,:])
Y.append(data[1,:])

X = np.column_stack(X)
Y = np.column_stack(Y)

# plot data
plt.plot( X, Y, 'b.')
plt.plot( X,-Y, 'b.')
plt.plot(-X, Y, 'b.')
plt.plot(-X,-Y, 'b.')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
