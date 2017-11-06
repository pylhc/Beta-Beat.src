from __future__ import print_function
import __init__
import os, sys
import numpy as np
import pandas as pd
from math import factorial
from Utilities import tfs_pandas


################################
#           CLASSES
################################

class TwissOptics(object):
    ''' 
    Generate all Hamiltonian terms contributions 
    per element for AC dipole excitations 
    '''

    def __init__(self, model_path):
        self.mad_twiss = tfs_pandas.read_tfs(model_path)
        self.plane = 'H'
        self.phase_advances = self._get_phase_advances()

    def _get_phase_advances(self):
        ''' Calculate phase advances between all elements ''' 
        self.mad_twiss.set_index('NAME', inplace=True)
        plane_mu = "MUX" if self.plane == "H" else "MUY"

        phases_mdl = np.array(self.mad_twiss.loc[self.mad_twiss.index, plane_mu])
        phase_advances = pd.DataFrame((phases_mdl[np.newaxis,:] - phases_mdl[:,np.newaxis]) % 1.0)
        phase_advances.set_index(self.mad_twiss.index, inplace=True)
        phase_advances.columns = self.mad_twiss.index
        return phase_advances

    @property
    def betas(self):
        ## BETAS FUNCTION
        return self.new_df.betas 

    @betas.setter
    def betas(self, new):
        self.new_df.betas = new 

    @property
    def dispersion(self):
        ## DISPERSION FUNCTION
        return self.new_df.dispersion

    @dispersion.setter
    def dispersion(self, new):
        self.new_df.dispersion = new 

    @property
    def coupling(self):
        ## COUPLING FUNCTION
        return self.new_df.coupling
    
    @coupling.setter
    def coupling(self, new):
        self.new_df.coupling = new 


class DrivingTerms(object, TwissOptics):

    def __init__(self, model_path):
        super().__init__(model_path)
        self._make_tune_step_mask()

    def fterm(self, orders):

        self.driving_terms = self._calc_terms_from_sources(self.order, )

    def _make_tune_step_mask(self):
        acd_idx = 100
        ''' Create a mask to destinguish between before and after AC dipole '''
        ref = np.zeros([len(self.mad_twiss), len(self.mad_twiss)])
        mask_A = np.tril(ref+1, -1)
        ref[acd_idx+1:, :acd_idx] = 1
        mask_B = ref
        mask_C = -np.transpose(mask_B)
        self.tune_step_mask = mask_A - mask_B + mask_C

    def _calc_terms_from_sources(order, machine):
        tw = self.mad_twiss
        Qx, Qy = tw.Qx, tw.Qy
        r = list(order)
        j, k, l, m = int(r[1]) , int(r[2]) , int(r[3]), int(r[4])  
        ## ADD PARSING OF ELEMENTS FOR SPECIFIC SOURCE
        factor = factorial(j) * factorial(k) * factorial(l) * factorial(m) * 2**(j+k+l+m)

        if (l+m)%2==0:
            src = 'K'+str(j+k+l+m-1)+'L'
            sources_df['h'+str(order)] = -1.0/factor * (1j**(l+m)).real * tw[src] * tw['BETX']**((j+k)/2.) * tw['BETY']**((l+m)/2.)
            sources_df['f'+str(order)] = sources_df['h'+str(order)] / (1. - np.exp(2 * np.pi * 1j*( (j-k)*Qx + (l-m)*Qy ) ))
        elif (l+m)%2==1:
            src = 'K'+str(j+k+l+m-1)+'SL'
            sources_df['h'+str(order)] = -1.0/factor * (1j**(l+m)).imag * tw[src] * tw['BETX']**((j+k)/2.) * tw['BETY']**((l+m)/2.) 
            sources_df['f'+str(order)] = sources_df['h'+str(order)] / (1. - np.exp(2 * np.pi * 1j*( (j-k)*Qx + (l-m)*Qy ) ))
        
        return sources_df 




################################
#           MAIN
################################

def do_main():
    ''' Main function launcher '''
    model_path = 'twiss_2octupoles.dat'
    TwissOptics(model_path)

    
if __name__ == '__main__':
    do_main()
