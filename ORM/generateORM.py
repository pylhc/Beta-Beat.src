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
def orbit(bpms, vc, twiss):	# calculate orbit response to dipole using analytical expression based on twiss parameters. Will meed to be adjusted for different machines, to determine which correctors are H or V
##############################

	ft=open('twiss.resp.dat','w')

	ft.write('* NAME S X Y QX QY\n')
	ft.write('$ %s %le %le %le %le %le\n')

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
		ft.write(bpm+' '+str(s)+' '+str(dx)+' '+str(dy)+' '+str(Qx)+' '+str(Qy)+'\n')
		u=u+1
	
	ft.close()

    
######################### For ORM
def writeparams(variable, increment):
#########################
    g = open ('changeparametersORM', 'w')
    g.write(vq+' = '+variable+' + '+str(increment)+';\n')
    g.close()
    return
    
########################
def justtwiss():
#########################
    global dictionary
    system('madx < job.ORM.madx > scum')
    x=twiss('twiss.orbit.dat')
    return x
    
################### Main

it=open('iteration.dat','r')
iteration=int(it.read())
	
FullResponse={}   #Initialize FullResponse

execfile('results/variableNames.py')
varQ=quadNames()          # magnet K1 error (thin quad in center of QF, QD)
varCH=corNamesH()
varCV=corNamesV()
varU=bpmNames()
print 'varQ = ', varQ
print 'varCH = ', varCH
print 'varCV = ', varCV
print 'varU = ', varU

execfile('results/AllCalib_'+str(iteration)+'.py')
calCH=corCalbH()
calCV=corCalbV()
calX=bpmCalbH()
calY=bpmCalbV()
print "dCalH = ",calCH
print "dCalV = ",calCV
print "xCal = ",calX
print "yCal = ",calY


incrQ=0.0001				# increment by which variable parameters will be changed
FullResponse['incrQ']=incrQ		#Store this info for future use
	

for vq in varQ: # loop over quads for strength; for dummy quad k0, loop over BPM calib

	writeparams(vq, incrQ)
	MADtwiss=justtwiss()

	d=0
	for vc in varCH:
		var=vq+'-'+vc
		orbit(varU,vc,MADtwiss)
		FullResponse[var]=twiss('twiss.resp.dat')
		print var
		d=d+1

	d=0
	for vc in varCV:
		var=vq+'-'+vc
		orbit(varU,vc,MADtwiss)
		FullResponse[var]=twiss('twiss.resp.dat')
		print var
		d=d+1
		   
pickle.dump(FullResponse,open('FullResponse','w'),-1)



