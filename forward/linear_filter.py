import numpy as np 
import matplotlib.pyplot as plt

class SLB():
    '''
    '''
    def __init__(self):
        self.__ab2 = None
        self.__rhoap_obs = None
        self.__rhoap_cal = None
        self.__rhotr = None
        self.__thick = None
        self.__rms_err = None
        self.__phi_r = None
        self.__a_r = None


    def run(self, ab2, rhoap_obs, rhotr, thick, filter_coeff='guptasarma_7'):
        '''
        parameter:
            :ab2 : half distance of electrode [m], 1D np.array
            :rhoap_obs : apperant resistivity from observation data [ohm.m], 1D np.array
            :rhotr : true resistivity of earth model, 1D np.array
            :thick : thicknes of layers, 1D np.array
            :filter_coeff : type of filter coefficien (e.g default 'guptasarma_7', sevent filter of guptasarma filter)
        
        '''
        self.__ab2 = ab2
        self.__rhoap_obs = rhoap_obs
        self.__rhotr = rhotr
        self.__thick = thick
        self.__filter_coefficent(filter_coeff)
        self.__rhoap_cal = np.array([self.__lin_fil(ab) for ab in self.__ab2])
        self.__rms_err = self.__rms_error()
        return self.__rhoap_cal

    def plot_mod(self):
        '''
            for plotting of curve matching and earth model
        parameter:
            save_fig : if 'True' will saving figure.
        '''
        d = []
        d.append(0)
        for i in range(len(self.__thick)):
            d.append(d[-1] + self.__thick[i])

        d = np.repeat(d, 2)

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        ax[0].set_title('Curve Matching')
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
        for i in range(len(self.__rhotr)-1):
            img = np.hstack((img, np.repeat(self.__rhotr[i], 100*self.__thick[i])))
        img = np.array(img[1:].reshape(-1, 1))

        L, R, T, B = min(self.__rhotr), max(self.__rhotr), min(d[1:]), max(d[1:])
        border = [ L+ 0.001*L, R+ 0.02*R, T+ 0.02*T, B+ 0.02*B]
        im = ax[1].imshow(img, aspect='auto', cmap='jet', origin='lower', extent=border)

        ax[1].set_title('Earth Model')
        ax[1].plot(np.repeat(self.__rhotr, 2)[:-1], d[1:], color='w', linewidth=2 )
        ax[1].set_xscale('log')
        ax[1].set_xlabel(r'$\rho_{true}  [\Omega m]$') 
        ax[1].set_ylabel('Depth [m]')
        ax[1].invert_yaxis()   
        fig.colorbar(im,ax=ax[1], label=r'$\rho_{true}  [\Omega m]$')
        return fig

    def __lin_fil(self, L):
        layer = len(self.__rhotr) - 1
        T = self.__rhotr[-1] 

        T_phi = []
        for i in range(len(self.__a_r)):
            lambda_r = 10 ** (self.__a_r[i] - np.log10(L))
            while layer > 0:
                layer -=1
                tanh = np.tanh(lambda_r*self.__thick[layer])
                T = (T + self.__rhotr[layer] * tanh) / (1 + T*tanh/self.__rhotr[layer])
        
            T_phi.append(T*self.__phi_r[i])

        rho_app = np.sum(T_phi)
        return rho_app

    def __rms_error(self):
        err = np.sqrt(np.mean((self.__rhoap_obs - self.__rhoap_cal) ** (2)))
        return err

    def __filter_coefficent(self, filter_coeff):
        if filter_coeff == 'guptasarma_7':
            a = np.array([-0.17445, 0.09672, 0.36789, 0.63906, 0.91023, 1.1814, 1.45257])
            phi = ([0.1732, 0.2945, 2.147, -2.1733, 0.6646, -0.1215, 0.0155])
        
        elif filter_coeff == 'guptasarma_11':
            a = np.array([-0.420625, -0.20265625, 0.0153125, 0.23328125, 0.45125, 0.66921875, 0.8871875, 1.10515625, 1.323125, 1.54109375, 1.7590625])
            phi = np.array([0.041873, -0.022258, 0.38766, 0.647103, 1.84873, - 2.96084, 1.358412, -0.37759, 0.097107, -0.024243, 0.004046])
 
        elif filter_coeff == 'guptasarma_22':
            a = np.array([-0.980685, -0.771995, -0.563305, -0.354615, -0.145925, 0.062765, 0.271455, 0.480145,  0.688835,  0.897525,  1.106215, 1.314905, 1.523595, 1.732285, 1.940975, 2.149665,  2.358355,  2.567045, 2.775735])
            phi = np.array([0.00097112, -0.00102152, 0.00906965, 0.01404316, 0.09012,  0.30171582, 0.99627084, 1.3690832, -2.99681171, 1.65463068, -0.59399277,  0.22329813 , -0.10119309,  0.05186135, -0.02748647, 0.01384932, -0.00599074, 0.00190463, -0.0003216 ])
        
        self.__phi_r = phi
        self.__a_r = a

    @property
    def rms_err(self):
        return self.__rms_err

if __name__ == "__main__":
    fname = './samples/sample1.txt'
    data = np.loadtxt(fname)
    ab2 = data[:, 0]
    rhoap_obs = data[:, 1]
    rhotr = np.array([120, 30, 2])
    thick = np.array([10, 10])
    
    lf = SLB()
    lf.run(ab2, rhoap_obs, rhotr, thick, filter_coeff='guptasarma_7')
    lf.plot_mod()
    plt.show()
