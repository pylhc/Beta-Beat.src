#!/usr/bin/env pythonafs

from Numeric import *
from LinearAlgebra import *
import sys
from os import system
import math
from Numeric import *
#from pymadtable import madtable
from metaclass import twiss
from SLSknobs import *


#########################
def writeparams(deltafamilie):
#########################
    global variables
    g = open ('changeparameters', 'w')
    i=0
    for var in variables:
        g.write(var+' ='+ var+'+ ('+str(deltafamilie[i])+');\n')
        i +=1
    g.close()
    return



#########################
def justtwiss(deltafamilies):
#########################
    global dictionary
    print deltafamilies
    writeparams(deltafamilies)
    system('madx < sls_f6cwo.RM.madx > scum')
    x=twiss('twiss.dat',dictionary)
    return x





execfile('mydictionary.py')

FullResponse={}   #Initialize FullResponse
variables=varSLS()            #Define variables
delta1=zeros(len(variables))*1.0   #Zero^th of the variables
incr=ones(len(variables))*0.01    #increment of variables


FullResponse['incr']=incr           #Store this info for future use
FullResponse['delta1']=delta1       #"     "     "


for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        FullResponse[variables[i]]=justtwiss(delta)

FullResponse['0']=justtwiss(delta1) #Response to Zero, base , nominal

pickle.dump(FullResponse,open('FullResponse.Numeric','w'),-1)



