import sys
import matplotlib.pyplot as plt 
sys.path.append('../')

import numpy as np
import fitting
from fitting import SLB_DLSInv

fname = '../samples/sample1.txt'
data = np.loadtxt(fname)
ab2 = data[:, 0]
rhoap_obs = data[:, 1]

rhotr = np.array([10, 11, 30, 35])
thick = np.array([5, 15, 18])

inv = SLB_DLSInv()
rho, thick = inv.fit(ab2, rhoap_obs, rhotr, thick, dumping=0.1, epsilon=0.0005 , err_min= 0.000001, filter_coeff='guptasarma_7')

print('rho model :', rho)
print('thickness :', thick)
# inv.plot_err()
plt.colorbar()
plt.show()