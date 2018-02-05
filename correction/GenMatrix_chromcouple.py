r"""
.. module: PythonClasses4MAD.GenMatrix_chromcouple

Created on ??

TODO: description


.. moduleauthor:: Unknown
"""

import os
import datetime
import time

from numpy.linalg import pinv as generalized_inverse
from numpy import dot as matrixmultiply
import numpy as np


def make_list(x, m, modelcut, errorcut, CorD):
    """Makes list in coupling correction """
    if x == []:
        return []

    result_names_list = []
    count = 0
    for i in range(len(x.NAME)):
        bn = x.NAME[i].upper()
        if bn in m.indx:
            i_x = x.indx[bn]
            i_m = m.indx[bn]
            if ((abs(x.Cf1001r[i_x]-m.Cf1001r[i_m]) < modelcut)
                        and (abs(x.Cf1001i[i_x]-m.Cf1001i[i_m]) < modelcut)
                        and (x.Cf1001iERR[i_x] < errorcut)
                        and (x.Cf1001rERR[i_x] < errorcut)):
                result_names_list.append(x.NAME[i])
        else:
            print "Not in Response:", bn
            count += 1
    
    if count > 0:
        print "Warning: ", count, "BPMs removed from data for not beeing in the model"
    
    return result_names_list


def write_params(deltafamilie, variables, app=0, path="./"):
    if (app == 0):
        mode ='w'
    if (app == 1):
        mode = 'a'
    timestamp = datetime.datetime.fromtimestamp(time.time())
    knobs_file = open(os.path.join(path, "changeparameters_chromcouple"), mode)
    tfs_file = open(os.path.join(path, "changeparameters_chromcouple.tfs"), mode)
    print >>tfs_file, "@", "APP", "%le", app
    print >>tfs_file, "@", "PATH", "%s", path
    print >>tfs_file, "@", "DATE", "%s", timestamp.ctime()
    print >>tfs_file, "*", "NAME", "DELTA"
    print >>tfs_file, "$", "%s", "%le"
    for i, var in enumerate(variables):
        knobs_file.write(var+' = '+ var+' + ( '+str(deltafamilie[i])+' );\n')
        tfs_file.write(var+'   '+str(deltafamilie[i])+'\n')
    knobs_file.close()
    tfs_file.close()


def correctcouple(a, chromcouple_input, cut=0.01, app=0, path="./"):
    R = np.transpose(chromcouple_input.sensitivity_matrix)
    vector = chromcouple_input.compute_vector(a)
    wg = chromcouple_input.wg
    
    len_couplelist = len(chromcouple_input.couplelist)
    weisvec = np.array(np.concatenate(
                                      [np.sqrt(wg[0])*np.ones(len_couplelist),
                                      np.sqrt(wg[1])*np.ones(len_couplelist),
                                      np.sqrt(wg[2])*np.ones(len_couplelist),
                                      np.sqrt(wg[3])*np.ones(len_couplelist)]
                                      )
                      )
    
    Rnew = np.transpose(np.transpose(R)*weisvec)
    
    delta = -matrixmultiply(generalized_inverse(Rnew,cut), (vector-chromcouple_input.zerovector)/chromcouple_input.normvector)
    
    write_params(delta, chromcouple_input.varslist, app,  path=path)
    
    return [delta, chromcouple_input.varslist]



class ChromCoupleInput:
    def __init__(self, varslist, couplelist=[], wg=[1,1,1,1,1]):
        self.varslist = varslist
        self.couplelist = couplelist
        self.wg = wg
        self.sensitivity_matrix = []
        self.zerovector = []


    def compute_vector(self,a):
        Cf1001r = []
        Cf1001i = []
        Cf1010r = []
        Cf1010i = []
        
        for bpm_name in self.couplelist:
            Cf1001r.append(a.Cf1001r[a.indx[bpm_name]])
            Cf1001i.append(a.Cf1001i[a.indx[bpm_name]])
            Cf1010r.append(a.Cf1010r[a.indx[bpm_name]])
            Cf1010i.append(a.Cf1010i[a.indx[bpm_name]])
        
        return np.array(np.concatenate([Cf1001r,Cf1001i,Cf1010r,Cf1010i]))


    def computeSensitivityMatrix(self,x):
        #global zerovector, normvector
        self.zerovector = self.compute_vector(x['0'])
        incr=x['incr'][0]  # BUG! need to read it from FullResponse!
        
        self.normvector = np.ones(4*len(self.couplelist))
        for var in self.varslist:
            vector=self.compute_vector(x[var])
            
            self.sensitivity_matrix.append((vector-self.zerovector)/self.normvector/incr)
        self.sensitivity_matrix=np.array(self.sensitivity_matrix)
    
        return self.sensitivity_matrix


