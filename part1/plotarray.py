import numpy as np 
import matplotlib.pyplot as plt 
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable 
plt.rcParams.update({'font.size': 12})

init = np.loadtxt("init0.dat", delimiter=",")
exact = np.loadtxt("analytical0.dat", delimiter=',')

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)

IM = ax.imshow(init, cmap='jet', vmin=0.0, vmax=1.0)
CS = ax.contour(init, 4, colors='white')
ax.clabel(CS, inline=1, color='white', fontsize=12)
fig.colorbar(IM, cax=cax, orientation='vertical' )
ax.set_xlabel(r"$NX_{1}$")
ax.set_ylabel(r"$NX_{2}$")
ax.set_title("Initial Condition")
plt.show()