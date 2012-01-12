# This program will calculate the coupling 

# Author: Glenn Vanbavinckhove

from Numeric import *
from metaclass import *
from math import *
import os, sys

print 'Welcome to Getcoupling.py ....'
print 'This program will calculate the coupling for the specified files'


#file=('/afs/cern.ch/user/g/gvanbavi/Coupling/ALLBPMs.coup')
file=('/afs/cern.ch/user/r/rtomas/w1/MAD/LHC/Beta-Beat/SimulateData/MCBCV.6R1.B1')

# Try using /afs/cern.ch/user/r/rtomas/w1/MAD/LHC/Beta-Beat/SimulateData/ALLBPMs.coup

fx=twiss(file+'_linx')
fy=twiss(file+'_liny')


# looking for x and y at the same plane
sx=fx.S
sy=fy.S


# Find double plane BPMs
print 'Finding dual BPMs'
DoubleBPMs=filter(lambda x: x in fx.NAME   , fy.NAME)
print 'Number of dual bpms:',len(DoubleBPMs)



# calculating the RDT (difference)
print ''
print 'Calulating the Resonance driving terms'
C10=fx.AMP01
C01=fy.AMP10
ftotal=[]    #array(len(C10),Int)


for bpm  in DoubleBPMs:

	f=0.5*atan(sqrt(C10[fx.indx[bpm]]*C01[fy.indx[bpm]]))

        #print 'BPM',bpm,'the value is',f
	
	print 'The RDT for ',bpm,'is', f

       	ftotal.append(f)

# calculating q

print '\n Start calculating q'

phi10=fx.MUX
theta10=fy.PHASE10
phi01=fx.PHASE01
theta01=fy.MUY


q1=[]
q2=[]

for bpm in DoubleBPMs:

	#print phi01[fx.indx[bpm]],'and', theta01[fy.indx[bpm]]
	
	q1=phi10[fx.indx[bpm]]-theta10[fy.indx[bpm]]-math.pi/2
	q2=phi01[fx.indx[bpm]]-theta01[fy.indx[bpm]]-math.pi/2
	
	print 'q1 is ', q1, 'q2 is', q2

	
# the global coupling coefficient

print '\n Start calculating the global coupling coefficient', '\n First calculating the amplitude'

Qx=fx.Q1  # Get tunes like this
Qy=fy.Q2

SUM=sum(ftotal)
CoupAmp=4*abs(Qx-Qy)*(SUM/len(DoubleBPMs))
print '\t => The amplitude is', CoupAmp,'\n'

print '\n =>Start calculating the phase'

phase=(1/len(DoubleBPMs))*sum(q1-phi10-theta01)+math.pi*1

print '\t => the phase is',phase










