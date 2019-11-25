import numpy as np 
import matplotlib.pyplot as plt 
import math 
plt.rcParams.update({'font.size': 12})

error = np.loadtxt("PJacobi2078.dat", delimiter=",")
plt.imshow(error)
plt.colorbar(label=r"$|T-T_{\rm sol}|$")
plt.xlabel(r"$NX_{1}$")
plt.ylabel(r"$NX_{2}$")
plt.title("Point-SOR-Jacobi Error")
plt.show()