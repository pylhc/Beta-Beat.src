# Created 28/02/2009
# This python script is a copy of the scripts for SPS and LHC
# Adjustment are made to fit SOLEIL
# Adjustments made by GLenn Vanbavinckhove


# General inputs
from Numeric import *
from LinearAlgebra import *
import sys
from os import system
import math
from Numeric import *
from metaclass import twiss
from AllLists import *

# Initialising variables
FullResponseQ={} # full response matrix quadrupole
FullResponseS={} # full response matrix sextupole
# => quadrupole variables
variablesQ = ["q1","q2","q3","q4","q5","q6","q7","q8","q9","q10"]
# => sextupole variables
variablesS = ["s1","s2","s3","s4","s5","s6","s7","s8","s9","s10"]

# => write parameters quadrupole and sextupole
#########################
def writeparamsQ(deltafamilie,magnet):
#########################
    global variables
    #if magnet=="S":
    g = open ('changeparameters'+'S', 'w')
    #elif magnet=='Q':
        #g = open ('changeparameters'+'Q', 'w')
    i=0
    for var in variables:
        g.write(var+' ='+ var+'+ ('+str(deltafamilie[i])+');\n')
        i +=1
    g.close()
    return

# => Perform madx job
#########################
def justtwiss(deltafamilies,magnet):
#########################
    global dictionary
    print deltafamilies
    writeparams(deltafamilies,magnet)
    if(magnet=="Q"):
        system('madx < job.iterateQ.madx > scum')
       
    elif(magnet=='S'):
        system('madx < job.iterateS.madx > scum')
               
    x=twiss('twiss_Response.dat')
    x.Cmatrix()
    return x

# = > Main method

   # variables for main method
bumpStrength=0.000088  # check this for SOLEIL
delta1=zeros(len(variables))
incr=ones(len(variables))*bumpStrength

   # Fullresponse matrix for quadrupoles
FullResponseQ['incr']=incr # storing the increment done
FullResponseQ['delta1']=delta1 # change applied
FullResponseQ['0']=justtwiss(delta1) # Response matrix to NO change in the variables (base)

   # Fullresponse matrix for sextupoles
FullResponseS['incr']=incr # storing the increment done
FullResponseS['delta1']=delta1 # change applied
FullResponseS['0']=justtwiss(delta1) # Response matrix to NO change in the variables (base)

# !! # loop implementing change in quadrupole magnet (variables)
   # => for quadrupoles
for i in range(0,len(delta1)) :
    delta=array(delta1)
    delta[i]=delta[i]+incr[i]
    FullResponseQ[variables[i]]=justtwiss(delta,Q)

   #=> dumping repsonse matrix for quadrupoles
pickle.dump(FullResponseQ,open('FullResponseQ.Numeric','w'),-1)

print 'Full response matrix for quadrupoles finished \n will proceed with sextupoles now ...'
   
# !! # loop implementing change in sextupole magnet (variables)
   # => for sextupoles
for i in range(0,len(delta1)) :
    delta=array(delta1)
    delta[i]=delta[i]+incr[i]
    FullResponseQ[variables[i]]=justtwiss(delta,Q)

   #=> dumping repsonse matrix for sextupoles
pickle.dump(FullResponseS,open('FullResponseS.Numeric','w'),-1)

print 'Full response matrix for sextupole finished \n Python job finished'
