import numpy as np 
import matplotlib.pyplot as plt 
import math
import pandas as pd 

xcoord = np.linspace(0,1.0,100)
pec = 20.0
y = (np.exp(xcoord*pec)-1.0)/(math.exp(pec)-1.0)

data = pd.read_csv("testFTBCS630.csv")

plt.plot(xcoord, y)
plt.plot(data.coord, data.u)
plt.show()