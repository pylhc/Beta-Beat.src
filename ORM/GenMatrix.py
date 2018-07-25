#!/usr/bin/env pythonafs


import sys
from os import system

from numpy import *
#from LinearAlgebra import *

import pickle

from metaclass import twiss
from variableNames import *
import datetime
import time


################################
def MakeList(measORM, calcORM, varslist):
################################
	if measORM==[]:
		return []

	t=[]
	cou=0
	print('measORM.name = ',measORM.NAME)
	for i in range(len(measORM.NAME)):
		sname=measORM.NAME[i].split('-')
		cor=sname[0]
		bpm=sname[1]
		val='kK0-'+cor
		try:
			calcORM[val]
		except:
			cou=cou+1
			print(val+' not present')
		else:
			m0=calcORM[val]
			try:
				m0.indx[bpm]
			except:
				cou=cou+1
				print(val+' not present')
			else:
				t.append(measORM.NAME[i])
	if cou > 0:
		print "Warning: ", cou, "correctors removed from data for not being in the model"

	return t			
	
############################################################
def writeparams(deltafamily, variables, oldCal, iteration, app, path):
##############################################################
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    
    a = datetime.datetime.fromtimestamp(time.time())
    itern=int(iteration)
    
    g = open ('results/madCalib.dat', mode)
    h = open ('results/AllCalib_'+str(itern+1)+'.py', mode)
    
    print 'len(deltafamily)=',len(deltafamily)
    print 'len(oldCal)=',len(oldCal)
     
    #print new calibration values to madCalib (quad only) and AllCalib (all):
    groupNames=['quadCalb()', 'corCalbH()', 'corCalbV()','bpmCalbH()','bpmCalbV()']
    
    j=0
    k=0
    print 'variables = ',variables
    for group in groupNames:
    	h.write('def '+group+':\n\tvar=[\n')
    	for var in variables[k]:
   	   	newcal=oldCal[j]+deltafamily[j]
   	 	#newcal=oldCal[j]*deltafamily[j]
   	 	h.write('\t'+str(newcal)+',\n')
		if k==0:
			g.write(variables[k][j]+' := '+str(newcal)+';\n')
		j=j+1
    	h.write('\t]\n\n\treturn var \n\n')
	    
    	k=k+1
    
   	    
    g.close()
    h.close()
    return

####################################################################
def correctbeatWei(ORMcalc, ORMmeas, minWei, beat_input, oldCal,  cut, iteration, app, path):
######################################################################
	
	wei=beat_input.weightsList(ORMmeas, minWei)
	print 'wei = ',wei
	
	R=   transpose(beat_input.response_matrix)
	print "dim R = ",shape(R)
	print "dim wei = ",shape(wei)
	RW=transpose(transpose(R)*wei)

	D=[]
	i=0
	for vec in transpose(RW):
		dLine=[]
		for j in range(len(transpose(RW))):
			if j==i:
				dLine.append(linalg.norm(vec))
			else:
				dLine.append(0.0)
		D.append(dLine)
		i=i+1
	print "dim D = ",shape(D)
	print "D = ",D
	
	vector=beat_input.computevectorEXP(ORMcalc, ORMmeas)
	vWei=vector*wei
	delta=dot(linalg.pinv(RW,cut),vWei) # sign - is removed for constructing error model
        
	'''
	if iteration==0:
		k=open(path+'sensM'+str(iteration)+'.dat','w')
   		i=0
		j=0
		for dim1 in R:
    		    for dim2 in R[i]:
		    	tmp=R[i]
			k.write(str(tmp[j])+"\t")
			j=j+1
		    i=i+1
		    j=0	
		    k.write("\n")
       		k.close()
       
	if iteration==0:
		k=open(path+'sensMW'+str(iteration)+'.dat','w')
   		i=0
		j=0
		for dim1 in RW:
    		    for dim2 in RW[i]:
		    	tmp=RW[i]
			k.write(str(tmp[j])+"\t")
			j=j+1
		    i=i+1
		    j=0	
		    k.write("\n")
       		k.close()
	
	
	k=open(path+'dVW'+str(iteration)+'.dat','w')
   	i=0
	for dim1 in vWei:
    	    k.write(str(vWei[i])+"\n")
	    i=i+1
	k.close()
	
	k=open(path+'dV'+str(iteration)+'.dat','w')
   	i=0
	for dim1 in vector:
    	    k.write(str(vector[i])+"\n")
	    i=i+1
	k.close()
		
	 
	vWnew=dot(RW,delta)-vWei
	k=open(path+'dVWpred_'+str(iteration+1)+'.dat','w')
   	i=0
	for dim1 in vWnew:
    	    k.write(str(vWnew[i])+"\n")
	    i=i+1
	k.close()
	'''
	
	writeparams(delta, beat_input.varslist, oldCal,  iteration, app, path)
	return [delta, beat_input.varslist]

#########################################
class beat_input:
##########################################
	def __init__(self, varslist, mlist=[], wg=[1]):
			self.varslist=varslist
			self.mlist=mlist
			self.wg=wg

	#############################################################################
	def computevectorEXP(self,ORMcalc,ORMmeas):	# make vector of differences between measured and model orbit response
	###############################################################################

    		orx=[] # Orbit response
		ory=[]

    		for m in self.mlist:
        		[cor,bpm]=m.split('-')
			val0='kK0-'+cor
			if cor[5]=='V':
				orx.append(0.0)
				ory.append(ORMmeas.MY[ORMmeas.indx[m]]-ORMcalc[val0].Y[ORMcalc[val0].indx[bpm]])
        		else:
				orx.append(ORMmeas.MX[ORMmeas.indx[m]]-ORMcalc[val0].X[ORMcalc[val0].indx[bpm]])
				ory.append(0.0)
				
		return array(concatenate([orx,ory]))
	
	#############################################################################
	def weightsList(self, ORMmeas, minWei):		# calculate weights for measured differences (based on uncertainty of measurement)
	###############################################################################

    		weiX=[] # uncertainty of measured Orbit response
		weiY=[]
				
		for m in self.mlist:
        		[cor,bpm]=m.split('-')
			if cor[5]=='V':
				if ORMmeas.dMY[ORMmeas.indx[m]]==0.0:
					weiX.append(0.0)
					weiY.append(0.0)
				else:
					weiX.append(0.0)
					weiY.append(1/max(ORMmeas.dMY[ORMmeas.indx[m]],minWei))
			else:
				if ORMmeas.dMX[ORMmeas.indx[m]]==0.0:
					weiX.append(0.0)
					weiY.append(0.0)
				else:
					weiX.append(1/max(ORMmeas.dMX[ORMmeas.indx[m]],minWei))
					weiY.append(0.0)
			print m, ', wei = ', weiX[-1],', dx = ', ORMmeas.dMX[ORMmeas.indx[m]]
		return array(concatenate([weiX,weiY]))
		
	#################################################################################
	def computevector(self,ORMcalc,knob,incr,group):
	##################################################################################
    		orx=[] # Orbit response
		ory=[]
		
		if group==0:  ## quad calibration
			quad=knob
			for m in self.mlist:
				[mcor,bpm]=m.split('-')
				val=quad+'-'+mcor
				val0='kK0-'+mcor
				try:
					ORMcalc[val]
				except:
					orx.append(0.0)
					ory.append(0.0)
					
				else:
					if mcor[5]=='V':
						orx.append(0.0)
						ory.append((ORMcalc[val].Y[ORMcalc[val].indx[bpm]]-ORMcalc[val0].Y[ORMcalc[val0].indx[bpm]])/incr)
					else:
						orx.append((ORMcalc[val].X[ORMcalc[val].indx[bpm]]-ORMcalc[val0].X[ORMcalc[val0].indx[bpm]])/incr)
						ory.append(0.0)
						
		if group==1:  ## H dipole calibration
			cor=knob
			for m in self.mlist:
				[mcor,bpm]=m.split('-')
				val0='kK0-'+mcor
				if cor==mcor:
					orx.append(-ORMcalc[val0].X[ORMcalc[val0].indx[bpm]])
					ory.append(0.0)
				else:
					orx.append(0.0)
					ory.append(0.0)
				
					
		if group==2:  ## V dipole calibration
			cor=knob
			for m in self.mlist:
				[mcor,bpm]=m.split('-')
				val0='kK0-'+mcor
				if cor==mcor:
					orx.append(0.0)
					ory.append(-ORMcalc[val0].Y[ORMcalc[val0].indx[bpm]])
				else:
					orx.append(0.0)
					ory.append(0.0)				
					
					
		if group==3:  ## H BPM calibration
			bpm=knob
			for m in self.mlist:
				[mcor,mbpm]=m.split('-')
				val0='kK0-'+mcor
				if mcor[5]=='V':
					orx.append(0.0)
					ory.append(0.0)
				elif bpm==mbpm:
					orx.append(ORMcalc[val0].X[ORMcalc[val0].indx[bpm]])
					ory.append(0.0)
				else:
					orx.append(0.0)
					ory.append(0.0)
			
		if group==4:  ## V BPM calibration
			bpm=knob
			for m in self.mlist:
				[mcor,mbpm]=m.split('-')
				val0='kK0-'+mcor
				if mcor[5]=='H':
					orx.append(0.0)
					ory.append(0.0)
				elif bpm==mbpm:
					orx.append(0.0)
					ory.append(ORMcalc[val0].Y[ORMcalc[val0].indx[bpm]])
				else:
					orx.append(0.0)
					ory.append(0.0)
			
		return array(concatenate([orx,ory]))


	###################################################################################
	def computeResponseMatrix(self,ORMcalc):
	################################################################################
    		self.response_matrix=[]
    		incrQ=ORMcalc['incrQ']
		print 'incrQ = ', incrQ
		
		#change in response from quad calibration:
		for quad in self.varslist[0]: 
			vector=self.computevector(ORMcalc, quad,incrQ,0)
			self.response_matrix.append(vector)
			
		# change in response from H dipole calibration:
		for dip in self.varslist[1]: 
			vector=self.computevector(ORMcalc, dip, 0, 1)
			self.response_matrix.append(vector)
			
		# change in response from V dipole calibration:
		for dip in self.varslist[2]: 
			vector=self.computevector(ORMcalc, dip, 0, 2)
			self.response_matrix.append(vector)
		
		# change in response from bpm calibration:
		for bpm in self.varslist[3]: 
			vector=self.computevector(ORMcalc, bpm, 0, 3)
			self.response_matrix.append(vector)
			
		# change in response from bpm calibration:
		for bpm in self.varslist[3]: 
			vector=self.computevector(ORMcalc, bpm, 0, 4)
			self.response_matrix.append(vector)
		
		self.response_matrix=array(self.response_matrix,dtype=float)
		return self.response_matrix

    

