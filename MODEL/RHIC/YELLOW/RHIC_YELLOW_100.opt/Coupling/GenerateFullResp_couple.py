#!/usr/bin/env pythonafs

from Numeric import *
from LinearAlgebra import *
import sys, math
from os import system
#from pymadtable import madtable
from metaclass import twiss
from AllLists_couple import *
import pickle

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
    system('madx < job.iterate_couple.madx > scum')
    x=twiss('twiss_couple.dat')
    print x.Q1, x.Q2
    return x






FullResponse={}   #Initialize FullResponse
#variables= pickle.load(open('qCirc.pic','r'))         #Define variables

variables= varsRHIC()

delta1=zeros(len(variables))*1.0   #Zero^th of the variables
incr=ones(len(variables))*0.0001    #increment of variables


FullResponse['incr']=incr           #Store this info for future use
FullResponse['delta1']=delta1       #"     "     "
FullResponse['0']=justtwiss(delta1) #Response to Zero, base , nominal


for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        FullResponse[variables[i]]=justtwiss(delta)

pickle.dump(FullResponse,open('FullResponse.Numeric_couple','w'),-1)



