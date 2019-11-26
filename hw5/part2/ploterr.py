import numpy as np 
import matplotlib.pyplot as plt 
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable 
import pandas as pd
plt.rcParams.update({'font.size': 12})

w08 = pd.read_csv("lsorw08_res.log")
w10 = pd.read_csv("lsorw10_res.log")
w15 = pd.read_csv("lsorw12_res.log")

err_w08_after = w08.res[2:len(w08.itr)-1].values
err_w08_before = w08.res[1:len(w08.itr)-2].values

err_w10_after = w10.res[2:len(w10.itr)-1].values
err_w10_before = w10.res[1:len(w10.itr)-2].values

err_w15_after = w15.res[2:len(w15.itr)-1].values
err_w15_before = w15.res[1:len(w15.itr)-2].values

plt.plot(w08.itr[1:len(w08.itr)-2].values, np.log10(err_w08_after/err_w08_before), label=r"$\omega=0.8$")
plt.plot(w10.itr[1:len(w10.itr)-2].values, np.log10(err_w10_after/err_w10_before), label=r"$\omega=1.0$")
plt.plot(w15.itr[1:len(w15.itr)-2].values, np.log10(err_w15_after/err_w15_before), label=r"$\omega=1.2$")
plt.ylabel(r"residual")
plt.legend(loc='lower right')

'''
plt.plot(w08.itr, w08.err, label=r"$\omega=0.8$")
plt.plot(w10.itr, w10.err, label=r"$\omega=1.0$")
plt.plot(w15.itr, w15.err, label=r"$\omega=1.2$")
plt.legend(loc='upper right')
'''
plt.xlim(0,)
plt.xlabel("iteration")
plt.title("LSOR-GS")
plt.show()
