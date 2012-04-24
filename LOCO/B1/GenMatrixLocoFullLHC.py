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
from AllLists import *
from LinearAlgebra import *
import datetime
import time




################################
def MakeList(x,m,varslist):
################################
	if x==[]:
		return []

	vq=varslist[0]
	vc=varslist[1]
	fval=vq[0]+'-'+vc[0]
	m0=m[fval]

	t=[]
	cou=0
	for i in range(len(x.NAME)):
		sname=x.NAME[i].split('-')
		cor=sname[0]
		bpm=sname[1]
		val=vq[0]+'-'+cor
		m0=m[val]
		try:
			m0.indx[bpm]
		except:
			print "Not in Response:", cor.upper()
			cou=cou+1
		else:
			t.append(x.NAME[i])
	if cou > 0:
		print "Warning: ", cou, "correctors removed from data for not beeing in the model"

	
	return t
			

######################################
def MakeListOne(x,m,varslist,single):
######################################
	if x==[]:
		return []

	vq=varslist[0]
	vc=varslist[1]
	fval=vq[0]+'-'+vc[0]
	m0=m[fval]

	t=[]
	cou=0
	for i in range(len(x.NAME)):
		sname=x.NAME[i].split('-')
		cor=sname[0]
		bpm=sname[1]
		val=vq[0]+'-'+cor
		if vc[single]==cor:
			m0=m[val]
		else:
			m0=0
		try:
			m0.indx[bpm]
		except:
			print "Not in Response:", cor.upper()
			cou=cou+1
		else:
			t.append(x.NAME[i])
	if cou > 0:
		print "Warning: ", cou, "correctors removed from data for not beeing in the model"

	
	return t
			




############################################################
def writeparams(deltafamilie, variables, app=0, path="./"):
##############################################################
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open (path+'changeparameters', mode)
    f = open (path+'changeparameters.QdK', mode)
    print deltafamilie
    i=0
    for var in variables[0]:
        g.write(var+' = '+ var+' + ( '+str(deltafamilie[i])+' );\n')
        i +=1
    for var in variables[1]:
        #if 'ch' in var:
	#	g.write(var+'->HKICK := '+ var+'->HKICK * ( '+str(1.0+deltafamilie[i])+' );\n')
	#else:
	#	g.write(var+'->VKICK := '+ var+'->VKICK * ( '+str(1.0+deltafamilie[i])+' );\n')
	g.write(var+'->KICK := '+ var+'->KICK * ( '+str(1.0+deltafamilie[i])+' );\n')
        i +=1
    for var in variables[2]:
        g.write('!H '+var+' := '+ var+' * ( '+str(1.0+deltafamilie[i])+' );\n')
        i +=1
    for var in variables[2]:
        g.write('!V '+var+' := '+ var+' * ( '+str(1.0+deltafamilie[i])+' );\n')
        i +=1
    '''
    for var in variables[3]:
        g.write(var[0]+' = '+ var[0]+' + ( '+str(deltafamilie[i])+' );\n')
        i +=1
    for var in variables[4]:
        g.write(var[0]+' = '+ var[0]+' + ( '+str(deltafamilie[i])+' );\n')
        i +=1
    for var in variables[2]:
        g.write('!PSI '+var+' := '+ var+' * ( '+str(deltafamilie[i])+' );\n')
        i +=1
    for var in variables[1]:
        g.write('!PSI '+var+' := '+ var+' * ( '+str(deltafamilie[i])+' );\n')
        i +=1
    '''
    g.close()
    

    i=0
    for var in variables[5]:
        f.write(var+' '+str(-deltafamilie[i])+'\n') # Sign - for correction
        i +=1
    f.close()
    return




###########################################################
def correctbeat(a,beat_input, cut=0.01, app=0, path="./"):
##########################################################
	R=   transpose(beat_input.sensitivity_matrix) 
	vector=beat_input.computevector(a)
	wg=beat_input.wg
	weisvec=array(sqrt(wg[0])*ones(len(beat_input.mlist)))
	Rnew=transpose(transpose(R)*weisvec)
	delta=-matrixmultiply(generalized_inverse(Rnew,cut),vector/beat_input.normvector)
	writeparams(delta, beat_input.varslist, app,  path=path)
	return


####################################################################
def correctbeatEXP(ores,mor,beat_input, cut=0.01, app=0, path="./"):
######################################################################
	R=   transpose(beat_input.sensitivity_matrix)
        vector=beat_input.computevectorEXP(ores,mor)
        wg=beat_input.wg
	print len(R),len(beat_input.mlist)
	#if single==0:
	#weisvec=array(sqrt(wg[0])*ones(2*len(beat_input.mlist)))
	#else:
	#weisvec=array(sqrt(wg[0])*ones(2*len(beat_input.mlist)+1))
	#Rnew=transpose(transpose(R)*weisvec)
	Rnew=transpose(transpose(R)*1.0) # forget weighting for the time being
	#print (vector-beat_input.zerovector) # Good line to debug big correction strs
	delta=matrixmultiply(generalized_inverse(Rnew,cut),vector) # sign - is removed for constructing error model
        writeparams(delta, beat_input.varslist, app, path)
	return [delta, beat_input.varslist]




#########################################
class beat_input:
##########################################
	def __init__(self, varslist, mlist=[], wg=[1]):
			self.varslist=varslist
			self.mlist=mlist
			self.wg=wg


	#############################################################################
	def computevectorEXP(self,ores,mor):
	###############################################################################

    		orx=[] # Orbit response
		ory=[]

		

		incr=ores['incr'][0]
		qc=ores['qc']

    		for m in self.mlist:
        		[cor,bpm]=m.split('-')
			val0='q0-'+cor
        		orx.append(mor.MX[mor.indx[m]]-ores[val0].X[ores[val0].indx[bpm]]/(incr/qc))
			ory.append(mor.MY[mor.indx[m]]-ores[val0].Y[ores[val0].indx[bpm]]/(incr/qc))
    		return array(concatenate([orx,ory]))

	#################################################################################
	def computevector(self,ores,knob,single):
	##################################################################################
    		orx=[] # Orbit response
		ory=[]

		incr=ores['incr'][0]
		qc=ores['qc']
		#qe=ores['qe']

		if single==0:
			kq=knob
			for m in self.mlist:
				[cor,bpm]=m.split('-')
				val=kq+'-'+cor
				val0='q0-'+cor
				orx.append((ores[val].X[ores[val].indx[bpm]]-ores[val0].X[ores[val0].indx[bpm]])/(incr/qc))
				ory.append((ores[val].Y[ores[val].indx[bpm]]-ores[val0].Y[ores[val0].indx[bpm]])/(incr/qc))
		if single==1:
			cor=knob
			for m in self.mlist:
				[mcor,bpm]=m.split('-')
				val0='q0-'+mcor
				if cor==mcor:
					orx.append(ores[val0].X[ores[val0].indx[bpm]]/(incr/qc))
					ory.append(ores[val0].Y[ores[val0].indx[bpm]]/(incr/qc))
				else:
					orx.append(0.0)
					ory.append(0.0)
		if single==2:
			bpm=knob
			for m in self.mlist:
				[cor,mbpm]=m.split('-')
				val0='q0-'+cor
				if bpm==mbpm:
					orx.append(ores[val0].X[ores[val0].indx[bpm]]/(incr/qc))
					ory.append(0.0)
				else:
					orx.append(0.0)
					ory.append(0.0)
		if single==3:
			bpm=knob
			for m in self.mlist:
				[cor,mbpm]=m.split('-')
				val0='q0-'+cor
				if bpm==mbpm:
					orx.append(0.0)
					ory.append(ores[val0].Y[ores[val0].indx[bpm]]/(incr/qc))
				else:
					orx.append(0.0)
					ory.append(0.0)
		if single==4:
			sext=knob
			for m in self.mlist:
				[cor,bpm]=m.split('-')
				val=sext+'-'+cor
				val0='q0-'+cor
				orx.append((ores[val].X[ores[val].indx[bpm]]-ores[val0].X[ores[val0].indx[bpm]])/(incr/qc))
				ory.append((ores[val].Y[ores[val].indx[bpm]]-ores[val0].Y[ores[val0].indx[bpm]])/(incr/qc))
		if single==5:
			quad=knob
			for m in self.mlist:
				[cor,bpm]=m.split('-')
				val=quad+'-'+cor
				val0='q0-'+cor
				orx.append((ores[val].X[ores[val].indx[bpm]]-ores[val0].X[ores[val0].indx[bpm]])/(incr/qc))
				ory.append((ores[val].Y[ores[val].indx[bpm]]-ores[val0].Y[ores[val0].indx[bpm]])/(incr/qc))
		if single==6:
			bpm=knob
			for m in self.mlist:
				[cor,mbpm]=m.split('-')
				val0='q0-'+cor
				if bpm==mbpm:
					orx.append(ores[val0].Y[ores[val0].indx[bpm]]/(incr/qc)) # The sign seems OK to be consistent with the definition in the machine
					ory.append(ores[val0].X[ores[val0].indx[bpm]]/(incr/qc))
				else:
					orx.append(0.0)
					ory.append(0.0)
		if single==7:
			cor=knob
			for m in self.mlist:
				[mcor,bpm]=m.split('-')
				val0='q0-'+mcor
				#if 'ch' in mcor:
				#	mcoro=mcor.replace('ch','cv')
				#elif 'cv' in mcor:
				#	mcoro=mcor.replace('cv','ch')
				#val0='q0-'+mcoro
				if cor==mcor:
					if 'ch' in mcor:
						orx.append(0.0)
						ory.append(-ores[val0].Y[ores[val0].indx[bpm]]/(incr/qc))
					elif 'cv' in mcor:
						orx.append(ores[val0].X[ores[val0].indx[bpm]]/(incr/qc))
						ory.append(0.0)
				else:
					orx.append(0.0)
					ory.append(0.0)
		if single==999:
			cor=knob
			for m in self.mlist:
				[mcor,bpm]=m.split('-')
				val0='q0-'+mcor
				if cor==mcor:
					orx.append(ores[val0].X[ores[val0].indx[bpm]]/(incr/qc))
					ory.append(ores[val0].Y[ores[val0].indx[bpm]]/(incr/qc))
    		return array(concatenate([orx,ory]))


	###################################################################################
	def computeSensitivityMatrix(self,ores):
	################################################################################
    		self.sensitivity_matrix=[]
    		incr=ores['incr'][0]  # BUG! need to read it from FullResponse!
		qc=ores['qc']
		#qe=ores['qe']

		for kq in self.varslist[0]:
			print kq
			vector=self.computevector(ores,kq,0)
			self.sensitivity_matrix.append(vector/incr)
		for cor in self.varslist[1]: # Corrector calibration
			vector=self.computevector(ores,cor,1)
			self.sensitivity_matrix.append(vector)
		for bpm in self.varslist[2]: # BPM horizontal calibration
			print bpm
			vector=self.computevector(ores,bpm,2)
			self.sensitivity_matrix.append(vector)
		for bpm in self.varslist[2]: # BPM vertical calibration
			vector=self.computevector(ores,bpm,3)
			self.sensitivity_matrix.append(vector)
		'''
		#for sext in self.varslist[3]: # Sextupole misalignment (vertical)
		#	vector=self.computevector(ores,sext[0],4)
		#	self.sensitivity_matrix.append(vector/incr) # Unit is mrad as ignoring the factor qe=1000
		#for quad in self.varslist[4]: # Quad misalignment (rotation)
		#	vector=self.computevector(ores,quad[0],5) # Unit is mm as ignoring the factor qe=1000
		#	self.sensitivity_matrix.append(vector/incr)
		for bpm in self.varslist[2]: # BPM tilt
			vector=self.computevector(ores,bpm,6)
			self.sensitivity_matrix.append(vector)
		for cor in self.varslist[1]: # Corrector tilt
			vector=self.computevector(ores,cor,7)
			self.sensitivity_matrix.append(vector)
		'''
		#print self.sensitivity_matrix
    		self.sensitivity_matrix=array(self.sensitivity_matrix)
		return self.sensitivity_matrix
    

