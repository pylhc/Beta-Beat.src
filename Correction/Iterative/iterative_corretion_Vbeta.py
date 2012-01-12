### Module for Iterative correction
### V1.0 created by Masa. Aiba 10/Mar/2009

from metaclass import *
from Numeric import *
from math import *
import cmath
import sys, pickle,os
from os import system
#import operator
from string import *
from AllLists import *
from optparse import OptionParser



def GetPhases(prex,prey,j):
    phasex={}
    phasey={}
    for i in range(0,len(prex.NAME)-1):
        bn1=upper(prex.NAME[i])
        bn2=upper(prex.NAME2[i])
        phx=(j.MUX[j.indx[bn2]]-j.MUX[j.indx[bn1]])%1.0
        phasex[bn1]=phx
    for i in range(0,len(prey.NAME)-1):
        bn1=upper(prey.NAME[i])
        bn2=upper(prey.NAME2[i])
        phy=(j.MUY[j.indx[bn2]]-j.MUY[j.indx[bn1]])%1.0
        phasey[bn1]=phy

    return [phasex,phasey]



def phasesubtract(prex,prey,nit,options):


    MADTwiss1=twiss("./twiss.corrected.dat")
    [phase1x,phase1y]=GetPhases(prex,prey,MADTwiss1)
    fmodel=options.rpath+'/MODEL/'+options.ACCEL+'/twiss.dat'
    MADTwiss0=twiss(fmodel)
    [phase0x,phase0y]=GetPhases(prex,prey,MADTwiss0)
    fpresentx=options.path+'/getphasex.0'+str(nit)+'.out'
    fpresenty=options.path+'/getphasey.0'+str(nit)+'.out'
    faftx=open(fpresentx,'w')
    faftx.write('@ Q1 %le '+str(MADTwiss1.Q1-64)+'\n')
    faftx.write('* NAME   NAME2  S      PHASE  STDPH  PHXMDL\n')
    faftx.write('$ %s     %s     %le    %le    %le    %le\n')
    fafty=open(fpresenty,'w')
    fafty.write('@ Q2 %le '+str(MADTwiss1.Q2-59)+'\n')
    fafty.write('* NAME   NAME2  S      PHASE  STDPH  PHYMDL\n')
    fafty.write('$ %s     %s     %le    %le    %le    %le\n')

    for i in prex.NAME:
        i2=prex.NAME2[prex.indx[i]]
        try:
            phasex=prex.PHASE[prex.indx[i]]-(phase1x[i]-phase0x[i])
            faftx.write('"'+i+'" '+'"'+i2+'" '+str(MADTwiss0.S[MADTwiss0.indx[i]])+' '+str(phasex)+' '+str(prex.STDPH[prex.indx[i]])+' '+str(phase0x[i])+'\n')
            print i
        except:
            0.0
    faftx.close()
    system('cp '+fpresentx+' '+options.path+'/getphasex.out')

    for i in prey.NAME:
        i2=prey.NAME2[prey.indx[i]]
        #print i, MADTwiss0.S[MADTwiss0.indx[i]], prex.POS1[prex.indx[i]]
        try:
            phasey=prey.PHASE[prey.indx[i]]-(phase1y[i]-phase0y[i])
            #print (phase1[i][1]-phase0[i][1])
            fafty.write('"'+i+'" '+'"'+i2+'" '+str(MADTwiss0.S[MADTwiss0.indx[i]])+' '+str(phasey)+' '+str(prey.STDPH[prey.indx[i]])+' '+str(phase0y[i])+'\n')
        except:
            0.0

    fafty.close()
    system('cp '+fpresenty+' '+options.path+'/getphasey.out')





def runcorrection(nit,options):

    # Not clever implementation obliged by the fact that execfile("script.py") seems not to accept options for script.py.
    optlist=[['-a ',options.ACCEL],['-t ',options.TECH],['-n ',options.ncorr],['-p ',options.path],['-c ',options.cut],['-e ',options.errorcut],['-m ',options.modelcut],['-r ',options.rpath],['-s ',options.MinStr],['-j ',options.JustOneBeam]]
    runcorrect='pythonafs '+options.rpath+'/Correction/correct.py'
    for i in range(0,len(optlist)):
        runcorrect=runcorrect+' '+optlist[i][0]+str(optlist[i][1])
    system(runcorrect)

    
    fpknobs=options.path+'/changeparameters_all'
    preknobs=open(fpknobs,'r')
    predelta={}
    for line in preknobs:
        if line=='return;':
            break
        lines=line.split()
        kq=lines[0]
        predelta[str(kq)]=lines[5]
    preknobs.close()

    fknobs=options.path+'/changeparameters'
    knobs=open(fknobs,'r')
    delta={}
    for line in knobs:
        lines=line.split()
        kq=lines[0]
        delta[str(kq)]=lines[5]
    knobs.close()

    var=quadvarsb2()

    fpknobs=options.path+'/changeparameters_all'
    rmpknobs='rm '+fpknobs
    system(rmpknobs)


    newknobs=open(fpknobs,'w')
    for i in range(0,len(var)):
        try:
            predel=float(predelta[var[i]])
        except:
            predel=0.0
        try:
            curdel=float(delta[var[i]])
        except:
            curdel=0.0
        delall=predel-curdel # inverse correction
        newknobs.write(var[i]+' = '+var[i]+' + ( '+str(delall)+' );\n')
    newknobs.write('return;')
    newknobs.close()

    system('cp '+fpknobs+' ./')
    system('cp '+file1+' '+file1+'_'+str(nit))





########### START ###############


## Opions are basically the same to that of correct.py
## and are passed to correct.py.

parser = OptionParser()
parser.add_option("-a", "--accel", 
		 help="What accelerator: LHCB1 LHCB2 SPS RHIC",
		 metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-t", "--tech", 
		 help="Which algorithm: SVD MICADO",
		 metavar="TECH", default="SVD",dest="TECH")
parser.add_option("-n", "--ncorr", 
		 help="Number of Correctors for MICADO",
		 metavar="NCORR", default=5,dest="ncorr")
parser.add_option("-p", "--path",
		 help="Path to experimental files",
		 metavar="PATH", default="./",dest="path")
parser.add_option("-c", "--cut",
                  help="Singular value cut for the generalized inverse",
                  metavar="CUT", default=0.1 , dest="cut")
parser.add_option("-e", "--errorcut",
                  help="Maximum error allowed for the phase measurement",
                  metavar="ERRORCUT", default=0.013 , dest="errorcut")
parser.add_option("-m", "--modelcut",
                  help="Maximum difference allowed between model and measured phase",
                  metavar="MODELCUT", default=0.05 , dest="modelcut")
parser.add_option("-r", "--rpath",
                  help="Path to BetaBeat repository (default is the afs repository)",
                  metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction (default is 1e-6)",
                  metavar="MinStr", default=0.000001 , dest="MinStr")
parser.add_option("-j", "--JustOneBeam",
                  help="0 Just quads from one beam are used, 1 all quads (default is 0)",
                  metavar="JUSTONEBEAM", default=0 , dest="JustOneBeam")

###  Additional option for the iterative method

parser.add_option("-i", "--Iteration",
                  help="Number of iteration (default is 5)",
                  metavar="ITERATION", default=0 , dest="Iteration")
parser.add_option("-b", "--BigError",
                  help="Big errors (known by SBS approach) in MAD format file, e.g., kqtl7.r3b2 = kqtl7.r3b2 + ( +0.004 ); Note that the format is fixed, i. e., the fifth culumn should be dk.",
                  metavar="BigError", default="./" , dest="BE")


(options, args) = parser.parse_args()



# Save getphasex/y.out, getDx.out for correction

command='cp '+options.path+'/getphasex.out '+options.path+'/getphasex.00.out'
print command
try:
    system(command)
except:
    print 'No getphasex.out file in ',options.path
    print 'Correction is not possible... leaving Iterative Correction.'
    sys.exit()



command='cp '+options.path+'/getphasey.out '+options.path+'/getphasey.00.out'
try:
    system(command)
except:
    print 'No getphasey.out file in ',options.path
    print 'Correction is not possible... leaving Iterative Correction.'
    sys.exit()



command='cp '+options.path+'/getNDx.out '+options.path+'/getNDx.00.out'
fragndx=0
try:
    fndx=options.path+'/getNDx.out '
    ftemp=open(fndx,'r')
    system(command)
    fndx.close()
    fragndx=1
except:
    print 'No getNDx.out file in ',options.path
    print 'We continue correction ignoring the dispersion.'



# Initialize knobs with big errors if given (by SBS approach)
file1=options.path+'/changeparameters_all'
iniknobs=open(file1,'w')
if options.BE!="./":
    fbe=open(options.BE,'r')
    for line in fbe:
        iniknobs.write(line)
    fbe.close()
iniknobs.write('return;')
iniknobs.close()
system('cp '+file1+' ./')
system('cp '+file1+' ./changeparameters_0')
system('cp '+file1+' '+file1+'_0')


# Load measurement data
filex=options.path+'/getphasex.00.out'
prex=twiss(filex)
filey=options.path+'/getphasey.00.out'
prey=twiss(filey)

# Copy MAD script to the current directory temporarily since MAD output are generated in the currect directory.
# Is there any other cleverer way??
cpMAD='cp '+options.rpath+'/Correction/Iterative/job.iterative.correction.madx ./'
system(cpMAD)

# Line to run MAD
runMAD='madx < ./job.iterative.correction.madx > out.mou'
#runMAD='madx < ./job.iterative.correction.madx'




# Iteration loop
nlimit=int(options.Iteration) # No. of iterations

for nit in range(1,nlimit+1):
    system(runMAD)
    phasesubtract(prex,prey,nit,options)
    runcorrection(nit,options)
    cpfile='cp ./twiss.corrected.dat '+options.path+'/twiss.'+str(nit)+'.dat'
    system(cpfile)

# Output twiss.corrected.dat at the end of iteration
system(runMAD)
mvfile='mv ./twiss.corrected.dat '+options.path+'/twiss.final.dat'
system(mvfile)

# Remove MAD script
rmfile='rm ./job.iterative.correction.madx'
system(rmfile)
#system('rm ./out.mou')

# Restore getphasex.out and getphasey.out
command='cp '+options.path+'/getphasex.00.out '+options.path+'/getphasex.out'
system(command)
command='cp '+options.path+'/getphasey.00.out '+options.path+'/getphasey.out'
system(command)
