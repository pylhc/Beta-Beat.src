#!/usr/bin/env pythonafs

import sys
from os.path import abspath, join, dirname, pardir
sys.path.append(abspath(join(dirname(__file__), pardir)))
from os import system

from numpy import *
from numpy.linalg import pinv as  generalized_inverse
from numpy import dot as matrixmultiply
import pickle

from Python_Classes4MAD.metaclass import twiss

##############################
def orbit(bpms, vc, twiss,twissall):
##############################
    
	u=0	
	for bpm in bpms:
		Qx=twiss.Q1
		Qy=twiss.Q2
		if vc[5]=='V':	# this condition will need to be changed if the corrector has a naming convention different from that in PSB
			dy=(1+calY[u])/(1+calCV[d])*sqrt(twiss.BETY[twiss.indx[bpm]]*twiss.BETY[twiss.indx[vc]])/(2*sin(pi*Qy))*cos(2*pi*abs(twiss.MUY[twiss.indx[bpm]]-twiss.MUY[twiss.indx[vc]])-pi*Qy)
			dx=0
		else:
			dx=(1+calX[u])/(1+calCH[d])*sqrt(twiss.BETX[twiss.indx[bpm]]*twiss.BETX[twiss.indx[vc]])/(2*sin(pi*Qx))*cos(2*pi*abs(twiss.MUX[twiss.indx[bpm]]-twiss.MUX[twiss.indx[vc]])-pi*Qx)
			dy=0.0
		s=twiss.S[twiss.indx[bpm]]
		ft.write(vc+'-'+bpm+' '+str(s)+' '+str(dx)+' '+str(dy)+' '+str(Qx)+' '+str(Qy)+'\n')
		u=u+1	
	     
	for j in twissall.NAME:
		Qx=twiss.Q1
		Qy=twiss.Q2
		if j.split('_')[0]=='DRIFT':
			continue
		if vc[5]=='V':	# this condition will need to be changed if the corrector has a naming convention different from that in PSB
			dy=sqrt(twissall.BETY[twissall.indx[j]]*twissall.BETY[twissall.indx[vc]])/(2*sin(pi*Qy))*cos(2*pi*abs(twissall.MUY[twissall.indx[j]]-twissall.MUY[twissall.indx[vc]])-pi*Qy)
			dx=0
		else:
			dx=sqrt(twissall.BETX[twissall.indx[j]]*twissall.BETX[twissall.indx[vc]])/(2*sin(pi*Qx))*cos(2*pi*abs(twissall.MUX[twissall.indx[j]]-twissall.MUX[twissall.indx[vc]])-pi*Qx)
			dy=0.0
		s=twissall.S[twissall.indx[j]]
		label=twissall.NAME[twissall.indx[j]]
		gt.write(vc+'-'+label+' '+str(s)+' '+str(dx)+' '+str(dy)+' '+str(Qx)+' '+str(Qy)+'\n')
	    
 

#########################
def justtwiss():
#########################
    global dictionary
    system('madx < job.ORM.madx > scum')
    x=twiss('twiss.orbit.dat')
    y=twiss('twiss.all.dat')
    return x, y


################### Main

it=open('iteration.dat','r')
iteration=int(it.read())

#read in the lists of BPM and corrector device names:
execfile('results/variableNames.py')
varCH=corNamesH()
varCV=corNamesV()
varU=bpmNames()
#print 'varQ = ', varQ
#print 'varCH = ', varCH
#print 'varCV = ', varCV
#print 'varU = ', varU

execfile('results/AllCalib_'+str(iteration+1)+'.py')
calCH=corCalbH()
calCV=corCalbV()
calX=bpmCalbH()
calY=bpmCalbV()
#print "calCH = ",calCH
#print "calCV = ",calCV
#print "calX = ",calX
#print "calY = ",calY


ft=open('results/ORM_calc_'+str(iteration+1)+'.dat','w')
ft.write('* NAME S X Y QX QY\n')
ft.write('$ %s %le %le %le %le %le\n')

gt=open('results/ORM_calc_All_'+str(iteration+1)+'.dat','w')
gt.write('* NAME S X Y QX QY\n')
gt.write('$ %s %le %le %le %le %le\n')

g = open ('changeparametersORM', 'w')		# remove old values from this file
g.close()
MADtwiss=justtwiss()
    
d=0	    
for vc in varCH:
	orbit(varU,vc,MADtwiss[0],MADtwiss[1])
	d=d+1

d=0	    
for vc in varCV:
	orbit(varU,vc,MADtwiss[0],MADtwiss[1])
	d=d+1

ft.close()
gt.close()

###################


# copy twiss file for each iteration (if needed for debugging)
system('cp twiss.all.dat results/twissAll_'+str(iteration+1)+'.dat')
    
