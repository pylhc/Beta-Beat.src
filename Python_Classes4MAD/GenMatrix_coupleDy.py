r"""
.. module: PythonClasses4MAD.GenMatrix_coupleDy

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
    if x==[]:
        return []

    result_names_list = []
    count = 0
    for i in range(len(x.NAME)):
        bn = x.NAME[i].upper()
        if bn in m.indx:
            i_x = x.indx[bn]
            i_m = m.indx[bn]
            if CorD == 'C':
                if ((abs(x.F1001W[i_x]-m.F1001W[i_m]) < modelcut) and (x.FWSTD1[i_x] < errorcut)):
                    result_names_list.append(x.NAME[i])
            elif CorD == 'D':
                if ((abs(x.DY[i_x]-m.DY[i_m]) < modelcut) and (x.STDDY[i_x] < errorcut)):
                    result_names_list.append(x.NAME[i])
        else:
            print "Not in Response:", bn
            count += 1
            
    if count > 0:
        print "Warning: ", count, "BPMs removed from data for not beeing in the model"
        
    return result_names_list


def write_params(deltafamilie, variables, app=0, path="./"):
    if (app == 0):
        mode = 'w'
    if (app == 1):
        mode = 'a'
    timestamp = datetime.datetime.fromtimestamp(time.time())
    knobs_file = open(os.path.join(path, "changeparameters_couple"), mode)
    tfs_file = open(os.path.join(path, "changeparameters_couple.tfs"), mode)
    print >>tfs_file, "@", "APP", "%le", app
    print >>tfs_file, "@", "PATH","%s", path
    print >>tfs_file, "@", "DATE", "%s", timestamp.ctime()
    print >>tfs_file, "*", "NAME", "DELTA"
    print >>tfs_file, "$", "%s", "%le"
    for i, var in enumerate(variables):
        knobs_file.write(var+' = '+ var+' + ( '+str(deltafamilie[i])+' );\n')
        tfs_file.write(var+'   '+str(deltafamilie[i])+'\n')
    knobs_file.close()
    tfs_file.close()


def coupling(a, b):
    # Applicable to both simulation-model and exp-model    
    rmsf1001 = np.sqrt(sum((a.F1001-b.F1001)**2)/len(b.F1001))
    rmsf1010 = np.sqrt(sum((a.F1010-b.F1010)**2)/len(b.F1001))
    peakf1001 = max(abs(a.F1001-b.F1001))
    peakf1010 = max(abs(a.F1010-b.F1010))
    return np.array([rmsf1001, rmsf1010, peakf1001, peakf1010])


def correctcouple(a, dispy, couple_input, cut=0.01, app=0, path="./"):
    R=   np.transpose(couple_input.sensitivity_matrix)
    vector=couple_input.computevector(a,dispy)
    wg=couple_input.wg
    
    weisvec = np.array(np.concatenate([
                                      np.sqrt(wg[0])*np.ones(len(couple_input.couplelist)),
                                      np.sqrt(wg[1])*np.ones(len(couple_input.couplelist)),
                                      np.sqrt(wg[2])*np.ones(len(couple_input.couplelist)),
                                      np.sqrt(wg[3])*np.ones(len(couple_input.couplelist)),
                                      np.sqrt(wg[4])*np.ones(len(couple_input.dispylist))
                                      ])
                      )
    Rnew = np.transpose(np.transpose(R)*weisvec)
    
    delta = -matrixmultiply(generalized_inverse(Rnew,cut), (vector-couple_input.zerovector)/couple_input.normvector)
    
    write_params(delta, couple_input.varslist, app,  path=path)
    
    return [delta, couple_input.varslist]



class CoupleInput:
    def __init__(self, varslist, couplelist=[], dispylist=[], wg=[1,1,1,1,1]):
        self.varslist = varslist
        self.couplelist = couplelist
        self.dispylist = dispylist
        self.wg = wg
        self.sensitivity_matrix = []
        self.zerovector = []

    def computevector(self,a,dispy):
        f1001r = []
        f1001i = []
        f1010r = []
        f1010i = []
        dy = []
        for bpm_name in self.couplelist:
            f1001r.append(a.F1001R[a.indx[bpm_name]])
            f1001i.append(a.F1001I[a.indx[bpm_name]])
            f1010r.append(a.F1010R[a.indx[bpm_name]])
            f1010i.append(a.F1010I[a.indx[bpm_name]])
        for bpm_name in self.dispylist:
            dy.append(dispy.DY[dispy.indx[bpm_name]])
        
        return np.array(np.concatenate([f1001r,f1001i,f1010r,f1010i,dy]))

    def computeSensitivityMatrix(self,x):
        self.zerovector = self.computevector(x['0'],x['0'])
        incr=x['incr'][0]  # BUG! need to read it from FullResponse!
        
        #ncouple=4*len(self.couplelist)
        self.normvector = np.array(np.concatenate([np.ones(len(self.couplelist)),np.ones(len(self.couplelist)),np.ones(len(self.couplelist)),np.ones(len(self.couplelist)),np.ones(len(self.dispylist)) ]))*1.0
        for var in self.varslist:
            if var in x:
                vector=self.computevector(x[var], x[var])
                self.sensitivity_matrix.append((vector-self.zerovector)/self.normvector/incr)
            else:
                raise KeyError("Variable "+var+" is not in Fullresponse.")
                
        self.sensitivity_matrix = np.array(self.sensitivity_matrix)
        
        return self.sensitivity_matrix


