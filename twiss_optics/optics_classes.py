from __future__ import print_function
import __init__
import os, sys
import numpy as np
import pandas as pd
from math import factorial
from Utilities import tfs_pandas
import matplotlib.pyplot as plt

################################
#           CLASSES
################################

class TwissOptics(object):

    def __init__(self, model_path):
        self.mad_twiss = tfs_pandas.read_tfs(model_path)
        self.mad_twiss.set_index('NAME', inplace=True)
        self.plane = 'H'
        self._phase_advances_x = self._get_phase_advances()
        self.plane = 'V'
        self._phase_advances_y = self._get_phase_advances()
        self._calc_terms_from_sources('f2002')
        self.f2002 = self.mad_twiss['f2002']
        
        self.f2002.abs().plot()
        plt.show()

    def _get_phase_advances(self):
        ''' Calculate phase advances between all elements ''' 
        plane_mu = "MUX" if self.plane == "H" else "MUY"
        
        phases_mdl = np.array(self.mad_twiss.loc[self.mad_twiss.index, plane_mu])
        phase_advances = pd.DataFrame((phases_mdl[np.newaxis,:] - phases_mdl[:,np.newaxis]) % 1.0)
        phase_advances.set_index(self.mad_twiss.index, inplace=True)
        phase_advances.columns = self.mad_twiss.index
        return phase_advances

    @property
    def betas(self):
        ## BETAS FUNCTION
        return self._results_df.betas 

    @betas.setter
    def betas(self, new):
        self._results_df.betas = new 

    @property
    def dispersion(self):
        ## DISPERSION FUNCTION
        return self._results_df.dispersion

    @dispersion.setter
    def dispersion(self, new):
        self._results_df.dispersion = new 

    @property
    def coupling(self):
        ## COUPLING FUNCTION
        return self._results_df.coupling
    
    @coupling.setter
    def coupling(self, new):
        self._results_df.coupling = new 

    def _calc_terms_from_sources(self, order):
        ## ADD PARSING OF ELEMENTS FOR SPECIFIC SOURCE AT START FOR tw = ....
        tw = self.mad_twiss
        Qx, Qy = tw.Q1, tw.Q2
        r = list(order)
        j, k, l, m = int(r[1]) , int(r[2]) , int(r[3]), int(r[4])  
        factor = factorial(j) * factorial(k) * factorial(l) * factorial(m) * 2**(j+k+l+m)
        phx = self._phase_advances_x.applymap(lambda x: np.exp(2 * np.pi * 1j * (j-k)*x))
        phy = self._phase_advances_y.applymap(lambda y: np.exp(2 * np.pi * 1j * (l-m)*y))
        phase_term = phx*phy

        if (l+m)%2==0:
            src = 'K'+str(j+k+l+m-1)+'L'
            beta_term = -1.0/factor * (1j**(l+m)).real * tw[src] * tw['BETX']**((j+k)/2.) * tw['BETY']**((l+m)/2.)
        elif (l+m)%2==1:
            src = 'K'+str(j+k+l+m-1)+'SL'
            beta_term = -1.0/factor * (1j**(l+m)).imag * tw[src] * tw['BETX']**((j+k)/2.) * tw['BETY']**((l+m)/2.) 
       
        local_hterms = phase_term.mul(beta_term)
        self.mad_twiss[order] = local_hterms.sum(axis=1)
        self.mad_twiss[order] = self.mad_twiss[order] / (1. - np.exp(2 * np.pi * 1j*( (j-k)*Qx + (l-m)*Qy ) ))
        


# class DrivingTerms(object, TwissOptics):
# 
#     def __init__(self, model_path):
#         super().__init__(model_path)
#         self._make_tune_step_mask()
# 
#     def fterm(self, orders):
# 
#         self.driving_terms = self._calc_terms_from_sources(self.order, )
# 
#     def _make_tune_step_mask(self):
#         acd_idx = 100
#         ''' Create a mask to destinguish between before and after AC dipole '''
#         ref = np.zeros([len(self.mad_twiss), len(self.mad_twiss)])
#         mask_A = np.tril(ref+1, -1)
#         ref[acd_idx+1:, :acd_idx] = 1
#         mask_B = ref
#         mask_C = -np.transpose(mask_B)
#         self.tune_step_mask = mask_A - mask_B + mask_C



################################
#           MAIN
################################

def do_main():
    ''' Main function launcher '''
    model_path = 'twiss_2octupoles.dat'
    two_octupoles = TwissOptics(model_path)


if __name__ == '__main__':
    do_main()
