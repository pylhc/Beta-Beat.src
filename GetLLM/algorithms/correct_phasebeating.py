'''
Created on May 5, 2017

@author: awegsche
'''
import re
import numpy as np

def get_correction_from_phasebeating(phase, commonbpms, twiss_elements, getllm_d):
      #---- create list of Ks, i.e. assign to each BPM a vector with all the errore lements that come after the bpm
    # and their respective errors
    # and model phases so that list_of_Ks[n] yields the error elements between BPM[n] and BPM[n+1]
    # update 2016-07-28: list_of_Ks[n][k], n: BPM number, k=0: quadrupole field errors,
    # k=1: transversal sextupole missalignments
    # k=2: longitudinal quadrupole missalignments  
    
    list_of_Ks = []
    
    for n in range(len(commonbpms) + getllm_d.range_of_bpms + 1):
        index_n = twiss_elements.indx[commonbpms[n % len(commonbpms)][1]]
        index_nplus1 = twiss_elements.indx[commonbpms[(n + 1) % len(commonbpms)][1]]
              
        quads = []
        magnetre = re.compile("^MQ")
               
        if index_n < index_nplus1:
            for i in range(index_n + 1, index_nplus1):
                if magnetre.match(twiss_elements.NAME[i]):
                    quads.append(i)
        list_of_Ks.append[quads]
        
    width = getllm_d.range_of_bpms / 2
    left_bpm = range(-width, 0)
    right_bpm = range(0 + 1, width + 1)
        
    BBA_combo = [[x, y] for x in left_bpm for y in left_bpm if x < y]
    ABB_combo = [[x, y] for x in right_bpm for y in right_bpm if x < y]
    BAB_combo = [[x, y] for x in left_bpm for y in right_bpm]
     
    combos = [[x,y] for x in range(getllm_d.range_of_bpms) for y in range(getllm_d.range_of_bpms) if x<y]

    
    
            
    A = np.zeros((len(combos), leng))
    print ":.............: dimensions of A:", A.shape
    
    for k in len(combos):
        pos = 0
        for j in xrange(combos[k][0]):
            pos += len(list_of_Ks[j][0])
        for j in xrange(combos[k][0], combos[k][1]):
            for l in xrange(len(list_of_Ks)):
                A[k][pos + j] = twiss_elements.BETX[list_of_Ks[j][0][l]]
            pos += len(list_of_Ks[j][0])
    
        