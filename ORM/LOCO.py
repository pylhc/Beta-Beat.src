import sys
sys.path.append('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD/')

import pickle
import operator 
from numpy import *
from os import system
from metaclass import twiss
import random,re,sys
#from LinearAlgebra import *
from optparse import OptionParser
from GenMatrix import *


########### START ###############

parser = OptionParser()

parser.add_option("-f", "--file",
		 help="Name of measured ORM data file",
		 metavar="FILE", default="measuredORM.dat",dest="file")
parser.add_option("-p", "--path",
		 help="Path to experimental data files",
		 metavar="PATH", default="./",dest="path")
parser.add_option("-c", "--cut",
                  help="Singular value cutoff for the pseudo-inverse",
                  metavar="SVCUT", default=0.0005 , dest="svcut")
parser.add_option("-e", "--errorcut",
                  help="Minimum uncertainty allowed for the orbit measurement",
                  metavar="ERRORCUT", default="0.001" , dest="errorcut")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction (default is 0.0001)",
                  metavar="MinStr", default=0.0001 , dest="MinStr")

(options, args) = parser.parse_args()
MinStr = float(options.MinStr)
minWei=float(options.errorcut) 
svcut= float(options.svcut)

datafilename = options.path+options.file
print 'ORM file = ', datafilename
ORMmeas=twiss(datafilename)

system('rm -r results')
system('mkdir results')
system('cp madCalib_0.dat results/madCalib.dat')
system('cp AllCalib_0.py results/AllCalib_0.py')
system('cp variableNames.py results/variableNames.py')

g = open ('iteration.dat', 'w')
g.write(str(-1))
g.close()   

#generate the model orbit response matrix, calculated analytically from MADX twiss output:
execfile('calcORM.py')

# one could change the number of iterations of the LOCO minimization procedure here:
for iteration in range(0,4):
    
    print 'iteration #',iteration
    g = open ('iteration.dat', 'w')
    g.write(str(iteration))
    g.close()
    
    print "Generating initial orbit response..."
    execfile('generateORM.py')
    FullResponse=pickle.load(open('FullResponse','r'))
    
    calbFile='results/AllCalib_'+str(iteration)+'.py'
    execfile(calbFile)
    
    calList=[quadCalb(),corCalbH(), corCalbV(),bpmCalbH(),bpmCalbV()]
    oldCalb=reduce(operator.add,calList)
    #print "oldCal= ",oldCalb  
    #print "len = ",len(oldCalb)   
     
    varslist=[quadNames(),corNamesH(),corNamesV(),bpmNames(),bpmNames()]
    #print "varslist= ",varslist
    #print "len varslist= ",len(varslist[0])," + ",len(varslist[1]) ," + ",len(varslist[2]) ," + ",len(varslist[3])," + ",len(varslist[4])   
	
    mlist=MakeList(ORMmeas, FullResponse, varslist)
    print "mlist = ", mlist
    #print 'len mlist = ',len(mlist)
    
    beat_inp=beat_input(varslist, mlist)
    sensitivity_matrix=beat_inp.computeResponseMatrix(FullResponse)
    [deltas, varslist] = correctbeatWei(FullResponse, ORMmeas, minWei, beat_inp, oldCalb, svcut, iteration, 0, "results/")
    print deltas

    execfile('calcORM.py')
   
