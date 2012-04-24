
import sys
sys.path.append('/tools/slsbd/TBTsrc/Python_Classes4MAD/')
sys.path.append('/tools/slsbd/TBTsrc/Numeric-23_p2.3/lib/python2.3/site-packages/Numeric/')

import pickle
from Numeric import *
from os import system
from metaclass import twiss
import random,re,sys
from LinearAlgebra import *
from optparse import OptionParser
from GenMatrixLocoFullLHC import *


########### START ###############


parser = OptionParser()
parser.add_option("-a", "--accel", 
		 help="What accelerator: SLS etc.",
		 metavar="ACCEL", default="SLS",dest="ACCEL")
parser.add_option("-p", "--path",
		 help="Path to experimental files",
		 metavar="PATH", default="./",dest="path")
parser.add_option("-f", "--file",
		 help="Orbit response file name",
		 metavar="DATAF", default="orbit.dat",dest="DATAF")
parser.add_option("-c", "--cut",
                  help="Singular value cut for the generalized inverse",
                  metavar="CUT", default=0.02 , dest="cut")
#parser.add_option("-e", "--errorcut",
#                  help="Maximum error allowed for the orbit measurement",
#                  metavar="ERRORCUT", default="0.013,0.2" , dest="errorcut")
parser.add_option("-m", "--modelcut",
                  help="Maximum difference allowed between model and measured phase",
                  metavar="MODELCUT", default="5.0" , dest="modelcut")
parser.add_option("-r", "--rpath",
                  help="Path to LOCO repository (default is the afs repository)",
                  metavar="RPATH", default="/tools/slsbd/LOCO" , dest="rpath")
parser.add_option("-o", "--optics",
                  help="optics",
                  metavar="OPT", default="/" , dest="OPT")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction (default is 0.0001)",
                  metavar="MinStr", default=0.0001 , dest="MinStr")
parser.add_option("-w", "--weight",
                help="Weighting factor",
                metavar="WGT", default="1", dest="WGT")


(options, args) = parser.parse_args()

    

print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
locopath = options.rpath
print "Path to Repository:", locopath
accelpath=locopath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/'
print "Path to Accelerator model", accelpath

datafile = options.DATAF

MinStr = float(options.MinStr)


##### implementing cuts for
    
modelcut=float(options.modelcut.split(",")[0])
#errorcut=float(options.errorcut.split(",")[0])
cut= float(options.cut)
print "Model and SVD cuts:", modelcut, cut
print "Starting loading Full Response optics"
FullResponse=pickle.load(open(accelpath+'/FullResponse.Numeric.MQ','r'))
FullResponseMQM=pickle.load(open(accelpath+'/FullResponse.Numeric.MQM','r'))
FullResponseMQT=pickle.load(open(accelpath+'/FullResponse.Numeric.MQT','r'))
FullResponseMQTL=pickle.load(open(accelpath+'/FullResponse.Numeric.MQTL','r'))
FullResponseMQTW=pickle.load(open(accelpath+'/FullResponse.Numeric.MQTW','r'))
FullResponseMQY=pickle.load(open(accelpath+'/FullResponse.Numeric.MQY','r'))
print "Loading ended"


## Merging the response matrices
keys=FullResponseMQM.keys()
for k in keys:
    FullResponse[k]=FullResponseMQM[k]

keys=FullResponseMQT.keys()
for k in keys:
    FullResponse[k]=FullResponseMQT[k]

keys=FullResponseMQTL.keys()
for k in keys:
    FullResponse[k]=FullResponseMQTL[k]

keys=FullResponseMQTW.keys()
for k in keys:
    FullResponse[k]=FullResponseMQTW[k]

keys=FullResponseMQY.keys()
for k in keys:
    FullResponse[k]=FullResponseMQY[k]


mor=twiss(options.path+'/'+datafile)
execfile(accelpath+'/AllLists.py')


varQ=MQ()+MQM()+MQT()+MQTL()+MQW()+MQY()
varCor=CorPartial()
varBPM=BPM()
varslist=[varQ,varCor,varBPM]




#for i in range(0,len(varC)):
wei=[]
for i in range(0,1):
    wei.append(int(options.WGT.split(",")[i]))

print "weight=",wei




mlist=MakeList(mor, FullResponse, varslist)
beat_inp=beat_input(varslist, mlist, wei)

print "Input ready"


single=0
sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)

TECH='SVD'
if TECH=="SVD": # Assume any other technique may be implemented in the future
    [deltas, varslist ] = correctbeatEXP(FullResponse, mor, beat_inp, cut=cut, app=0, path=options.path)
    ''' Skip iteration since it takes time. Later the computation time should be improved by improving the code and then turn it on.
    if 1:                           #All accelerators
        iteration=0
        # Let's remove too low useless correctors
        while (len(filter(lambda x: abs(x)< MinStr, deltas))>0):
            iteration=1+iteration
            il=len(varslist[0])
            varslist_t=[]
            for i in range(0,il):
              #print MinStr, deltas[i]
              if (abs(deltas[i]) > MinStr):
                varslist_t.append(varslist[0][i])
            varslist=[varslist_t,varslist[1]]
            if len(varslist[0])==0:
                print "You want to correct with too high cut on the corrector strength"
                sys.exit()
            beat_inp=beat_input(varslist, mlist, wei)
            sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
            [deltas, varslist] = correctbeatEXP(FullResponse, mor, beat_inp, cut=cut, app=0, path=options.path)
            print "Initial correctors:", il, ". Current: ",len(varslist), ". Removed for being lower than:", MinStr, "Iteration:", iteration
    '''
    print deltas
    
    

    
