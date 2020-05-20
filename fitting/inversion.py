import sys
sys.path.append('')

import numpy as np
import matplotlib.pyplot as plt

from forward import SLB
from numpy.linalg import inv, svd
from copy import deepcopy as dc

class SLB_LSInv():
    def __init__(self):
        self.__ab2 = None
        self.__rhoap_obs = None
        self.__rhoap_cal = None
        self.__rhotr_init = None
        self.__thick_init = None
        self.__epsilon = None
        self.__damping = None

        self.__rms_err = np.inf
        self.__err_hist = []
        self.__num_fig = None

        self.__M = None
        self.__N = None
        self.__m = None
        self.__n = None

        self.__slb = SLB()

    def fit(self, ab2, rhoap_obs, rhotr_init, thick_init, damping=0.01, epsilon=0.05, err_min=0.05, method='lm', filter_coeff='guptasarma_7'):

        #Check input parameters
        dim_err_mess = 'Dimension of {0} is invalid. Make sure {0} has 1D numpy array!'
        if ab2.ndim !=1:
            raise Exception(dim_err_mess.format('ab2')) 
        elif rhoap_obs.ndim !=1 :
            raise Exception(dim_err_mess.format('rhoap_obs'))
        elif rhotr_init.ndim != 1:
            raise Exception(dim_err_mess.format('rhotr_init'))
        elif thick_init.ndim !=1:
            raise Exception(dim_err_mess.format('thick_init'))

        if len(rhotr_init)-1 != len(thick_init):
            raise Exception('Your model is invalid. Number of thickness must be equel to (number of resistivities - 1)')
        
        if method != 'lm' and method != 'svd': 
            raise Exception('Method {0} is not supported! Please choose one `svd` or `lm` method!'.format(method))

        list_filter_coeff = ['guptasarma_7','guptasarma_11','guptasarma_22']
        if filter_coeff not in list_filter_coeff:
            raise Exception('Filter Coefficient {0} is not supported! Suported filter coefficient are {1}'.format(filter_coeff, list_filter_coeff))

        if epsilon <=0 or epsilon >=1:
            raise Exception('Epsilon range invalid! Epsilon must be greater than 0 and less than 1') 

        #Assign value of variable
        self.__ab2 = ab2
        self.__rhoap_obs = rhoap_obs
        self.__rhotr_init = rhotr_init
        self.__thick_init = thick_init
        self.__filter_coeff = filter_coeff
        self.__epsilon = epsilon
        self.__damping = damping
        self.__num_fig = 0

        self.__M = len(self.__ab2)          #number of data
        self.__n = len(self.__rhotr_init)   #number of layers
        self.__N = 2*self.__n - 1           #number of models

        self.__mod = np.vstack((self.__rhotr_init.reshape(-1, 1), self.__thick_init.reshape(-1, 1))) #models matricess
        self.__inversion(err_min, method)
        return self.__rhotr_init, self.__thick_init
           
    def __inversion(self, err_min, method):
        i = 0
        while self.__rms_err > err_min:
            self.__rhoap_cal = self.__slb.run(self.__ab2, self.__rhoap_obs, self.__rhotr_init, self.__thick_init, self.__filter_coeff)
            self.__rms_err = self.__slb.rms_err
            self.__err_hist.append( self.__rms_err)

            if i == 100 or self.__rms_err < err_min:
                break

            d = (self.__rhoap_obs - self.__rhoap_cal).reshape(-1, 1)          
            J = self.__jacobian()
            
            if method == 'lm':
                dmod = self.__lm_inv(J, d)
            
            elif method == 'svd':
                dmod = self.__svd_inv(J, d)

            self.__mod = self.__mod + dmod 

            self.__rhotr_init = self.__mod[:self.__n, 0]
            self.__thick_init = self.__mod[self.__n:, 0]
            i +=1

    def __lm_inv(self, J, d):
        Wm = self.__damping*np.eye(self.__N)
        dmod = inv(J.T @ J + Wm)@J.T@d
        return dmod

    def __svd_inv(self, J, d):
        U, S, Vh = svd(J, full_matrices=False)
        SS = S/(S + self.__damping)
        dmod = Vh.T@np.diag(SS)@U.T@d
        return dmod  


    def __jacobian(self):
        '''
        params:
            ab2 : half distance of electrode, 1D np.array
            rho : resistivities of layers, 1D np.array
            h   : thickness of layers, 1D np.array
            eps : partubation parameter, float
        return:
            J   : Jacobian matrices, 2D np.array
        '''
        del_m = self.__epsilon*(self.__mod.flatten())            #paturbation of models

        #calculate jacobian matrices using forward finite differences
        J = np.zeros((self.__M, self.__N))
        for j in range(self.__N):
            dm = dc(self.__mod.flatten())
            dm[j] = dm[j] + del_m[j]
            drho = dm[:self.__n]
            dh = dm[self.__n:]

            #forward finite differences
            f_mdm = self.__slb.run(self.__ab2, self.__rhoap_obs, drho, dh, self.__filter_coeff)
            f_m = self.__slb.run(self.__ab2, self.__rhoap_obs, self.__rhotr_init, self.__thick_init, self.__filter_coeff)
            J[:, j] = (f_mdm- f_m) / del_m[j]
        return J
    
    def plot_err(self):
        iter = np.arange(0, len(self.__err_hist))
        self.__num_fig += 1

        plt.figure(self.__num_fig)
        plt.plot(iter, self.__err_hist, marker='o')
        plt.ylabel('RMS Error [%]')
        plt.xlabel('Iterasi')
        plt.title('RMS Error')
    
    def plot_mod(self, save_fig = True):
        '''
            for plotting of curve matching and earth model
        parameter:
            save_fig : if 'True' will saving figure.
        '''
        self.__num_fig += 1

        d = []
        d.append(0)
        for i in range(len(self.__thick_init)):
            d.append(d[-1] + self.__thick_init[i])

        d = np.repeat(d, 2)

        _, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5), num = self.__num_fig)
        ax[0].set_title('Curve Fitting')
        ax[0].plot(self.__ab2, self.__rhoap_obs, color = 'r', marker='o', label= 'Data Observation ' )
        ax[0].plot(self.__ab2, self.__rhoap_cal, label= 'Data Forward')
        ax[0].set_xlabel(r'$\frac{AB}{2} [m]$')
        ax[0].set_ylabel(r'$\rho_{apperant}  [\Omega m]$')
        ax[0].set_yscale('log')
        ax[0].set_xscale('log')
        ax[0].grid(which='minor', color ='b', alpha=0.30)
        ax[0].grid(which='major', color ='r')
        ax[0].legend()

        img = [0]
        for i in range(len(self.__rhotr_init)-1):
            img = np.hstack((img, np.repeat(self.__rhotr_init[i], 100*self.__thick_init[i])))
        img = np.array(img[1:].reshape(-1, 1))

        L, R, T, B = min(self.__rhotr_init), max(self.__rhotr_init), min(d[1:]), max(d[1:])
        border = [ L+ 0.001*L, R+ 0.02*R, T+ 0.02*T, B+ 0.02*B]
        im = ax[1].imshow(img, aspect='auto', cmap='jet', origin='lower', extent=border)

        ax[1].set_title('Earth Model')
        ax[1].plot(np.repeat(self.__rhotr_init, 2)[:-1], d[1:], color='w', linewidth=2 )
        ax[1].set_xscale('log')
        ax[1].set_xlabel(r'$\rho_{true}  [\Omega m]$') 
        ax[1].set_ylabel('Depth [m]')
        ax[1].invert_yaxis()
        if save_fig : plt.savefig('Forward.png', dpi=120)
        plt.colorbar(im,ax=ax[1], label=r'$\rho_{true}  [\Omega m]$')


if __name__ == "__main__":
    fname = 'samples/sample1.txt'
    data = np.loadtxt(fname)
    ab2         = data[:, 0]
    rhoap_obs   = data[:, 1]

    rhotr = np.array([150, 20, 2])
    thick = np.array([20, 15])

    epsilon = 0.005
    err_min = 0.01
    damping = 0.01

    inversion = SLB_LSInv()
    rho, thick = inversion.fit(ab2, rhoap_obs, rhotr, thick, damping=damping, epsilon=epsilon, method='lm' , err_min= err_min, filter_coeff='guptasarma_22')

    print('rho model :', rho)
    print('thickness :', thick)
    
    inversion.plot_err()
    inversion.plot_mod()
    plt.show()
