#!/usr/bin/env pythonafs

from Numeric import *
from LinearAlgebra import *
import sys
from os import system
import math
from Numeric import *
from metaclass import twiss
#from AllLists import *


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
execfile('./AllLists.py')
exec('variables=Qb2()')           #Define variables
delta1=zeros(len(variables))*1.0   #Zero^th of the variables
incr=ones(len(variables))*0.0001    #increment of variables


FullResponse['incr']=incr           #Store this info for future use
FullResponse['delta1']=delta1       #"     "     "
#FullResponse['0']=justtwiss(delta1) #Response to Zero, base , nominal

#Set-up the madx file for

f=open('iter.madx','w')
for i in range(0,len(delta1)) : #Loop over variables
          delta=array(delta1)
          delta[i]=delta[i]+incr[i]
          var=variables[i]
          print >>f, var,"=", var, "+(",delta[i],");"
          print >>f, "twiss, file=twiss."+var+";"
#          print >>f, "system, \"gzip twiss."+var+" \";"
          print >>f, var,"=", var, "-(",delta[i],");"

print >>f, "twiss, file=twiss.0 ;"
f.close()

system('madx < job.iterate.madx')


for i in range(0,len(delta1)) : #Loop over variables
        delta=array(delta1)
        delta[i]=delta[i]+incr[i]
        var=variables[i]
        print "Reading twiss."+var
        FullResponse[var]=twiss("twiss."+var)
        #os.system('rm twiss.'+var)
FullResponse['0']=twiss('twiss.0') #Response to Zero, base , nominal

pickle.dump(FullResponse,open('FullResponse.Numeric','w'),-1)



