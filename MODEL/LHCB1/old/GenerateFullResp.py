#!/usr/bin/env pythonafs

from Numeric import *
from LinearAlgebra import *
import sys
from os import system
import math
from Numeric import *
#from pymadtable import madtable
from metaclass import twiss
from AllLists import *


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
    print deltafamilies
    writeparams(deltafamilies)
    system('madx < job.iterate.madx > scum')
    x=twiss('twiss.dat')
    return x






FullResponse={}   #Initialize FullResponse
variables=quadvarsb1()            #Define variables
delta1=zeros(len(variables))*1.0   #Zero^th of the variables
incr=ones(len(variables))*0.0001    #increment of variables


FullResponse['incr']=incr           #Store this info for future use
FullResponse['delta1']=delta1       #"     "     "
#FullResponse['0']=justtwiss(delta1) #Response to Zero, base , nominal


for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        FullResponse[variables[i]]=justtwiss(delta)

FullResponse['0']=justtwiss(delta1) #Response to Zero, base , nominal

pickle.dump(FullResponse,open('FullResponse.Numeric','w'),-1)



