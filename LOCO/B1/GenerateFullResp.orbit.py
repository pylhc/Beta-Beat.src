#!/usr/bin/env pythonafs

from Numeric import *
from LinearAlgebra import *
import sys
from os import system
import math
from Numeric import *
#from pymadtable import madtable
from metaclass import twiss


##############################
def orbit(bpms, cor, MADtwiss,ic,qc):
##############################

    ft=open('twiss.cod.dat','w')

    ft.write('* NAME S X Y\n')
    ft.write('$ %s %le %le %le\n')

    for bpm in bpms:
        if 'H' in cor:
            Q=MADtwiss.Q1
            bets=MADtwiss.BETX[MADtwiss.indx[bpm]]
            beti=MADtwiss.BETX[MADtwiss.indx[cor]]
            mus=MADtwiss.MUX[MADtwiss.indx[bpm]]
            mui=MADtwiss.MUX[MADtwiss.indx[cor]]
            kick=ic/qc
            dx=sqrt(bets)/2.0/sin(pi*Q)*kick*sqrt(beti)*cos(2.0*pi*abs(mus-mui)-pi*Q)
            dy=0.0
        else:
            Q=MADtwiss.Q2
            bets=MADtwiss.BETY[MADtwiss.indx[bpm]]
            beti=MADtwiss.BETY[MADtwiss.indx[cor]]
            mus=MADtwiss.MUY[MADtwiss.indx[bpm]]
            mui=MADtwiss.MUY[MADtwiss.indx[cor]]
            kick=ic/qc
            dy=sqrt(bets)/2.0/sin(pi*Q)*kick*sqrt(beti)*cos(2.0*pi*abs(mus-mui)-pi*Q)
            dx=0.0
        s=MADtwiss.S[MADtwiss.indx[bpm]]
        ft.write(bpm+' '+str(s)+' '+str(dx)+' '+str(dy)+'\n')

    ft.close()

    




######################### For LOCO
def writeparams(vq,incr):
#########################
    g = open ('changeparameters', 'w')
    g.write(vq+' = ' +vq+' + '+str(incr[0])+';\n')
    g.close()
    return

#########################
def justtwiss():
#########################
    global dictionary
    system('madx < job.twiss.madx > scum')
    #x=twiss('twiss.dat',dictionary)
    x=twiss('twiss.orbit.dat')
    return x



################### Main

#execfile('mydictionary.py')

twBPM=twiss('twiss.bpm.dat')
bpms=twBPM.NAME


execfile('AllLists.py')

FullResponse={}   #Initialize FullResponse
varQ=MQX()            #Define variables
varC=CorPartial()
delta1=zeros(len(varQ)*len(varC))*1.0   #Zero^th of the variables

qc=10
ic=0.0001
incr=ones(len(varQ)*len(varC))*ic    #increment of variables

FullResponse['qc']=qc
FullResponse['incr']=incr           #Store this info for future use
FullResponse['delta1']=delta1       #"     "     "


for vq in varQ : #Loop over variables
    writeparams(vq,incr)
    MADtwiss=justtwiss()
    for vc in varC:
        orbit(bpms, vc, MADtwiss,ic,qc)
        var=vq+'-'+vc
        FullResponse[var]=twiss('twiss.cod.dat')
        print var


pickle.dump(FullResponse,open('FullResponse.Numeric','w'),-1)



