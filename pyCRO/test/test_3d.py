'''
Description: 
Author: Hejun Xie
Date: 2020-10-11 10:56:01
LastEditors: Hejun Xie
LastEditTime: 2020-10-11 11:47:40
'''
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y, Z = axes3d.get_test_data(0.05)

# Plot the 3D surface
# ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)

print(Z.max())
print(Z.min())
# levels = np.linspace(-95, 100, 10)
levels = np.linspace(Z.min(), Z.max(), 100)
print(levels)

# Plot projections of the contours for each dimension.  By choosing offsets
# that match the appropriate axes limits, the projected contours will sit on
# the 'walls' of the graph
cset = ax.contourf(X, Y, Z, levels=levels, offset=0, cmap=cm.coolwarm)
# cset = ax.contourf(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
# cset = ax.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)

ax.set_xlim(-40, 40)
ax.set_ylim(-40, 40)
ax.set_zlim(-100, 100)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

