import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit

from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage.filters import gaussian_filter1d

import numpy as np
import math
#definir une fonction affine
def affine(x, a, b):
    return a*x+b;
x=[]
y=[]
z=[]
with open('/home/mat/tmp/point_sphere_uniform.csv', 'r') as csvfile:
    plots= csv.reader(csvfile, delimiter=';')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
        z.append(float(row[2]))


fig = plt.figure()


ax = fig.add_subplot((111),projection='3d')

ax.scatter(x, y, z, s=10,label='Uniforme 3D')
ax.set_aspect('equal')
ax.legend()
xyzlim = np.array([ax.get_xlim3d(),ax.get_ylim3d(),ax.get_zlim3d()]).T
XYZlim = [min(xyzlim[0]),max(xyzlim[1])]
ax.set_xlim3d(XYZlim)
ax.set_ylim3d(XYZlim)
ax.set_zlim3d(XYZlim)
plt.savefig('/home/mat/tmp/uniform3D_400pts.svg')

ax2d = fig.add_subplot((111))
ax2d.scatter(x, y,  s=10,label='Uniforme 3D, balayage')

ax2d.set_aspect('equal')
plt.savefig('/home/mat/tmp/uniform3D_400pts_balayage.svg')
plt.show()