#!/usr/bin/env pythonafs

from Numeric import *
from LinearAlgebra import *
import sys
from os import system
import math
from Numeric import *
#from pymadtable import madtable
from metaclass import twiss
#from AllLists import *
#from AllLists_couple import *



#########################
def writeparams(deltafamilie):
#########################
    global variables
    g = open ('changeparameters_couple', 'w')
    i=0
    for var in variables:
        g.write(var+' ='+ var+'+ ('+str(deltafamilie[i])+');\n')
        i +=1
    g.close()
    return



#########################
def justtwiss(deltafamilies):
#########################
    print deltafamilies
    writeparams(deltafamilies)
    system('madx < job.couplingb2.madx > scum')
    system('pythonafs twiss_converter_for_coupling.py twiss.couple.dat.org')
    x=twiss('twiss.couple.rm.dat')
    #x.Cmatrix()
    return x



##################
def squadvarsLHCB1():
    variables=[ 'kqs.r7b1 ',
                'kqs.r5b1 ',
                'kqs.r3b1 ',
                'kqs.r1b1 ',
                'kqs.l8b1 ',
                'kqs.l6b1 ',
                'kqs.l4b1 ',
                'kqs.l2b1 ',
                'kqs.a81b1 ',
                'kqs.a67b1 ',
                'kqs.a45b1 ',
                'kqs.a23b1 ']

    return variables


##################
def squadvarsLHCB2():
    variables=[ 'kqs.r2b2 ',
                'kqs.r4b2 ',
                'kqs.r6b2 ',
                'kqs.r8b2 ',
                'kqs.l1b2 ',
                'kqs.l3b2 ',
                'kqs.l5b2 ',
                'kqs.l7b2 ',
                'kqs.a12b2 ',
                'kqs.a78b2 ',
                'kqs.a56b2 ',
                'kqs.a34b2 ']

    return variables




FullResponse={}   #Initialize FullResponse
variables=squadvarsLHCB2()            #Define variables
delta1=zeros(len(variables))*1.0   #Zero^th of the variables
incr=ones(len(variables))*0.0001    #increment of variables


FullResponse['incr']=incr           #Store this info for future use
FullResponse['delta1']=delta1       #"     "     "
FullResponse['0']=justtwiss(delta1) #Response to Zero, base , nominal
system('mv twiss.couple.rm.dat twiss.dat') # The model twiss



f=open('resultsb2','w')

for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        x=justtwiss(delta)
        FullResponse[variables[i]]=x
        print >>f, variables[i], incr[i],x.F1001R[0], x.F1001I[0], x.F1010R[0], x.F1010I[0]

f.close()
#pickle.dump(FullResponse,open('FullResponse_couple.Numeric','w'),-1)

#system('rm twiss.couple.rm.dat')

