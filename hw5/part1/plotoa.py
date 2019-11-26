import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 12}) 

#pj = pd.read_csv("err_pj.csv")
#pj = pd.read_csv("err_pgs.csv")
#pj = pd.read_csv("err_lsor.csv")
pj = pd.read_csv("err_lsor.csv")

def straightline(x,x0,slope):
	return slope*(x-x0)

dx_ = pj.dx.iloc[2:]
dx_ = np.append(np.array([pj.dx[0]]), dx_)
l1_ = pj.L1[2:]
l1_ = np.append(np.array([pj.L1[0]]), l1_)



init_guess = [0,-2]
bestval1, covar_pj=curve_fit(straightline, np.log10(pj.dx[1:]), np.log10(pj.L1[1:]), p0=init_guess)
bestval2, covar_pj=curve_fit(straightline, np.log10(pj.dx[2:]), np.log10(pj.L1[2:]), p0=init_guess)
bestval3, covar_pj=curve_fit(straightline, np.log10(dx_), np.log10(l1_), p0=init_guess)
#bestval_pgs, covar_pgs=curve_fit(straightline, np.log10(pgs.dx[1:]), np.log10(pgs.L1[1:]), p0=init_guess)
print bestval1[1], bestval2[1], bestval3[1]

plt.scatter(np.log10(pj.dx.iloc[1:]), np.log10(pj.L1.iloc[1:]), c='r', label='__nolegend__')
#plt.scatter(np.log10(pgs.dx.iloc[1:]), np.log10(pgs.L1.iloc[1:]), c='b' )
#plt.scatter(np.log10(ljor.dx.iloc[1:]), np.log10(ljor.L1.iloc[1:]), c='g' )
#plt.scatter(np.log10(lsor.dx.iloc[1:]), np.log10(lsor.L1.iloc[1:]), c='k' )


xcoord = pj.dx.iloc[1:]
plt.plot(np.log10(xcoord), straightline(np.log10(xcoord), bestval1[0], bestval1[1]), 'k' ,label='method1, p='+str(bestval1[1]))
plt.plot(np.log10(xcoord), straightline(np.log10(xcoord), bestval2[0], bestval2[1]), 'k' ,label='method2, p='+str(bestval2[1]), ls='--')
plt.plot(np.log10(xcoord), straightline(np.log10(xcoord), bestval3[0], bestval3[1]), 'k' ,label='method3, p='+str(bestval3[1]), ls='-.')

plt.scatter(np.log10(pj.dx[0]), np.log10(pj.L1[0]), c='r', marker='+', label='__nolegend__')
#plt.scatter(np.log10(pgs.dx[0]), np.log10(pgs.L1[0]), c='b' , marker='+')
#plt.scatter(np.log10(ljor.dx[0]), np.log10(ljor.L1[0]), c='g', marker='+')
#plt.scatter(np.log10(lsor.dx[0]), np.log10(lsor.L1[0]), c='k' , marker='+')

plt.title(r"$o(\Delta)$"+ " LSOR-GS")
plt.legend(loc='upper left')
plt.xlabel(r"$\log_{10}\Delta$")
plt.ylabel(r"$\log_{10}L_{1}$")
plt.show()
