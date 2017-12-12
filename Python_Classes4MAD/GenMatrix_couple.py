#!/usr/bin/env pythonafs



# Just to make sure that the path to the libraires is defined 
import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy

import pickle
from Numeric import *
#from numpy.oldnumeric.linear_algebra import generalized_inverse
from os import system
from metaclass import twiss
import random,re
#from AllLists_couple import *
from LinearAlgebra import *
import datetime
import time
import cmath



################################
def MakeList(x,m):
################################
	t=[]
        cou=0
	if x==[]:
		return []
	for i in range(len(x.NAME)):
		try:
			m.indx[x.NAME[i].upper()]
		except:
			print "Not in Response:", x.NAME[i].upper()
			cou=cou+1
		else:
			t.append(x.NAME[i])
	if cou > 0:
		print "Warning: ", cou, "BPMs removed from data for not beeing in the model"
	return t	
			

########################
def MakePairs(x,m, modelcut=0.1, errorcut=0.027):
###########################################
	t=[]
        cou=0
	keys=x.__dict__.keys()
	if "PHYMDL" in keys:
		phmdl="PHYMDL"
	else:
		phmdl="PHXMDL"
	for i in range(len(x.NAME)):
		phm=x.__dict__[phmdl][i]
		
		if (x.STDPH[i] < errorcut and abs(x.PHASE[i]-phm) < modelcut):
			#print x.STDPH[i], errorcut, abs(x.PHASE[i]-phm)
			try:
				m.indx[x.NAME[i].upper()]
				m.indx[x.NAME2[i].upper()]
			except:
				print "Not in Response:", x.NAME[i].upper(), x.NAME2[i].upper()
				cou=cou+1
			else:
				t.append(x.NAME[i]+' '+x.NAME2[i])
			#print "Good:", x.NAME[i].upper(), x.NAME2[i].upper()
		else:
			cou=cou+1
	if cou > 0:
		print "Warning: ", cou, "BPM pairs removed from data for not beeing in the model or having too large error deviations: ", phmdl, modelcut, "STDPH",errorcut, "LEN", len(t)
	return t
		




##################################################
def writeparams(deltafamilie, variables, app=0, path="./"):
########################################################
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open (path+'changeparameters', mode)
    f = open (path+'changeparameters.tfs', mode)
    print >>f, "@", "APP", "%le", app
    print >>f, "@", "PATH","%s", path
    print >>f, "@", "DATE", "%s", a.ctime()
    print >>f, "*", "NAME", "DELTA"
    print >>f, "$", "%s", "%le"
    i=0
    print len(variables), len(deltafamilie)
    for var in variables:
        g.write(var+' = '+ var+' + ( '+str(deltafamilie[i])+' );\n')
	f.write(var+'   '+str(deltafamilie[i])+'\n')
        i +=1
    g.close();f.close()
    return




#################################################
def coupling(a,b):
#################################################
	# Applicable to both simulation-model and exp-model

	rmsf1001=sqrt(sum((a.F1001-b.F1001)**2)/len(b.F1001))
	rmsf1010=sqrt(sum((a.F1010-b.F1010)**2)/len(b.F1001))
	peakf1001=max(abs(a.F1001-b.F1001))
	peakf1010=max(abs(a.F1010-b.F1010))
	couplelist=MakePairs(a, b)
	return array([rmsf1001, rmsf1010, peakf1001, peakf1010])



##########################################
def correctcouple(a, dispy, couple_input, cut=0.01, app=0, path="./"):
################################################
	print "one"
	R=   transpose(couple_input.sensitivity_matrix) 
	vector=couple_input.computevector(a,dispy)
	wg=couple_input.wg
	print len(couple_input.couplelist), wg[2]
	weisvec=array(concatenate([sqrt(wg[0])*ones(len(couple_input.couplelist)),sqrt(wg[1])*ones(len(couple_input.couplelist)),sqrt(wg[2])*ones(len(couple_input.couplelist)),sqrt(wg[3])*ones(len(couple_input.couplelist)),sqrt(wg[4])*ones(len(couple_input.dispylist))]))
	Rnew=transpose(transpose(R)*weisvec)
	delta=-matrixmultiply(generalized_inverse(Rnew,cut),(vector-couple_input.zerovector)/couple_input.normvector)
	writeparams(delta, couple_input.varslist, app,  path=path)
	return [delta, couple_input.varslist]



##########################################
class couple_input:
##########################################
	def __init__(self, varslist, couplelist=[], dispylist=[], wg=[10,10,1,1,1]):
			self.varslist=varslist
			self.couplelist=couplelist
			self.dispylist=dispylist
			self.wg=wg


	##################################################
	def computevector(self,a,dispy):
	##################################################
		f1001r=[]
		f1001i=[]
		f1010r=[]
		f1010i=[]
	       	dy=[]
		for bpm in self.couplelist:
			f1001r.append(a.F1001R[a.indx[bpm]])
			f1001i.append(a.F1001I[a.indx[bpm]])
			f1010r.append(a.F1010R[a.indx[bpm]])
			f1010i.append(a.F1010I[a.indx[bpm]])
		for bpm in self.dispylist:
			dy.append(dispy.DY[dispy.indx[bpm]])

				
		return array(concatenate([f1001r,f1001i,f1010r,f1010i,dy]))



	###################################################################################
	def computeSensitivityMatrix(self,x):
	################################################################################
    		#global zerovector, normvector
    		self.zerovector=self.computevector(x['0'],x['0'])
    		self.sensitivity_matrix=[]
    		incr=x['incr'][0]  # BUG! need to read it from FullResponse!

    		#ncouple=4*len(self.couplelist)
    		self.normvector=array(concatenate([ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.dispylist)) ]))*1.0
    		for var in self.varslist:
        		vector=self.computevector(x[var], x[var])
        		self.sensitivity_matrix.append((vector-self.zerovector)/self.normvector/incr)
    		self.sensitivity_matrix=array(self.sensitivity_matrix)
		return self.sensitivity_matrix
    

########### START ###############
'''



print "Starting loading Full Response optics"
FullResponse=pickle.load(open('/afs/cern.ch/user/r/rtomas/w1/MAD/LHC/Beta-Beat/NewRespMatrixGen/FullResponse','r'))
print "Loading ended"

varslist=quadvarsb1()
variables=varslist
phasexlist=[]
phaseylist=[]
betaxlist=IndependentQuadsb1()
betaylist=betaxlist
displist=BPMsb1()
wei=[1,1,1,1,1,10] # Weights of phasex phasey betax betay disp and tunes

beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)



sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
    
y=twiss('twiss.base')

x=twiss('twiss.dat')

#--- evaluate beta-beat for output
print "Initial beatings:"
print "rmsx  rmsy  peakx  peaky  rmsphix  rmsphiy  peakphix  peakphiy   rmsdx  peakdx"
print betabeat(x,y)

#------ Apply svd correction
correctbeat(x, beat_inp, cut=0.04, app=0)
    
#----- compute twiss after correction
#system('madx < job.random.madx > scum  ')
#z=twiss('twiss.dat')
#bb1=betabeat(x,y);bb=betabeat(z,y)
    
'''
    
    

    
