### Module for Iterative correction
### V1.0 created by Masa. Aiba 15/Mar/2009
### V1.1 created by Glenn Vanbavinckhove 15/June/2009

from metaclass import *
from Numeric import *
from math import *
import cmath
import sys, pickle,os
from os import system
#import operator
from string import *
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



def phasesubtract(prex,prey,nit,options,path):

    global dictionary

    twissfile=path+"/twiss.corrected.dat"
    print "Where to read the twiss file "+twissfile
    MADTwiss1=twiss(twissfile,dictionary)
    [phase1x,phase1y]=GetPhases(prex,prey,MADTwiss1)
    fmodel=options.rpath+'/twiss.dat'
    MADTwiss0=twiss(fmodel,dictionary)
    [phase0x,phase0y]=GetPhases(prex,prey,MADTwiss0)
    fpresentx=options.path+'/getphasex.0'+str(nit)+'.out'
    fpresenty=options.path+'/getphasey.0'+str(nit)+'.out'
    faftx=open(fpresentx,'w')
    tunex=MADTwiss1.Q1 %1
    if tunex>0.5: tunex=tunex-1.0
    tuney=MADTwiss1.Q2 %1
    if tuney>0.5: tuney=tuney-1.0
    faftx.write('@ Q1 %le '+str(tunex)+'\n')
    faftx.write('* NAME   NAME2  S      PHASE  STDPH  PHXMDL\n')
    faftx.write('$ %s     %s     %le    %le    %le    %le\n')
    fafty=open(fpresenty,'w')
    fafty.write('@ Q2 %le '+str(tuney)+'\n')
    fafty.write('* NAME   NAME2  S      PHASE  STDPH  PHYMDL\n')
    fafty.write('$ %s     %s     %le    %le    %le    %le\n')

    for i in prex.NAME:
        i2=prex.NAME2[prex.indx[i]]
        try:
            phasex=prex.PHASE[prex.indx[i]]-(phase1x[i]-phase0x[i])
            faftx.write('"'+i+'" '+'"'+i2+'" '+str(MADTwiss0.S[MADTwiss0.indx[i]])+' '+str(phasex)+' '+str(prex.STDPH[prex.indx[i]])+' '+str(phase0x[i])+'\n')
        except:
            0.0
    faftx.close()
    system('cp '+fpresentx+' '+options.path+'/getphasex.out')

    for i in prey.NAME:
        i2=prey.NAME2[prey.indx[i]]
        #print i, MADTwiss0.S[MADTwiss0.indx[i]], prex.POS1[prex.indx[i]]
        try:
            phasey=prey.PHASE[prey.indx[i]]-(phase1y[i]-phase0y[i])
            print 'phase difference',(phase1y[i]-phase0y[i]) 
            #print (phase1[i][1]-phase0[i][1])
            fafty.write('"'+i+'" '+'"'+i2+'" '+str(MADTwiss0.S[MADTwiss0.indx[i]])+' '+str(phasey)+' '+str(prey.STDPH[prey.indx[i]])+' '+str(phase0y[i])+'\n')
        except:
            0.0

    fafty.close()
    system('cp '+fpresenty+' '+options.path+'/getphasey.out')





def runcorrection(nit,options,varlist):

    # Not clever implementation obliged by the fact that execfile("script.py") seems not to accept options for script.py.
    optlist=[['-a ',options.ACCEL],['-t ',options.TECH],['-n ',options.ncorr],['-p ',options.path],['-c ',options.cut],['-e ',options.errorcut],['-m ',options.modelcut],['-r ',options.rpath],['-s ',options.MinStr],['-j ',options.JustOneBeam]]
    runcorrect=options.pyth+' '+options.loc+'/Correction/correct.py'
    for i in range(0,len(optlist)):
        runcorrect=runcorrect+' '+optlist[i][0]+str(optlist[i][1])
    print runcorrect
    system(runcorrect)

    system('cp changeparameters changeparameters.save')

    
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


    fpknobs=options.path+'/changeparameters_all'
    rmpknobs='rm '+fpknobs
    system(rmpknobs)


    
    newknobs=open(fpknobs,'w')
    for i in range(0,len(varlist)):
        try:
            predel=float(predelta[varlist[i]])
        except:
            predel=0.0
        try:
            curdel=float(delta[varlist[i]])
        except:
            curdel=0.0
        delall=predel-curdel # inverse correction
        newknobs.write(varlist[i]+' = '+varlist[i]+' + ( '+str(delall)+' );\n')
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
                  metavar="MODELCUT", default=0.02 , dest="modelcut")
parser.add_option("-r", "--rpath",
                  help="Path to BetaBeat repository (default is the afs repository)",
                  metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/SPS/SPS" , dest="rpath")
parser.add_option("-l", "--loc",
                  help="location to beta beat",
                  metavar="LOC", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="loc")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction (default is 1e-6)",
                  metavar="MinStr", default=0.000001 , dest="MinStr")
parser.add_option("-j", "--JustOneBeam",
                  help="0 Just quads from one beam are used, 1 all quads (default is 0)",
                  metavar="JUSTONEBEAM", default=0 , dest="JustOneBeam")
parser.add_option("-x", "--python",
                  help="Which python to run",
                  metavar="PY", default="pytho" , dest="pyth")
parser.add_option("-v", "--Variables",
                  help="variables split with ,",
                  metavar="var", default="MQTb1" , dest="var")

###  Additional option for the iterative method

parser.add_option("-i", "--Iteration",
                  help="Number of iteration (default is 5)",
                  metavar="ITERATION", default=5 , dest="Iteration")
parser.add_option("-b", "--BigError",
                  help="Big errors (known by SBS approach) in MAD format file, e.g., kqtl7.r3b2 = kqtl7.r3b2 + ( +0.004 ); Note that the format is fixed, i. e., the sixth culumn should be dk.",
                  metavar="BigError", default="0" , dest="BE")
parser.add_option("-d", "--dictionary",
                help="File with the BPM dictionary",
                metavar="DICT", default="0", dest="dict")
parser.add_option("-z", "--mad",
                help="File with the BPM dictionary",
                metavar="MAD", default="/afs/cern.ch/group/si/slap/bin/", dest="mad")


(options, args) = parser.parse_args()


if options.Iteration=='0':
    print "Number of iteration cannot be equal to ZERO !! \n System will exit"
    sys.exit()


if options.dict!='0':
    execfile(options.dict)
else:
    dictionary={}

#### write directory file

filedir=options.path+"/pointer.madx"
writepoint=open(filedir,"w")
writepoint.write("system, 'ln -nfs "+options.rpath+" fileloc';\n")
writepoint.write("system, 'ln -nfs "+options.path+" filedes';\n")
writepoint.write("return;")
writepoint.close()


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
    fndx=options.path+'/getNDx.out'
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
if options.BE!="0":
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
# This avoids to output files in the repository.
# Is there any other cleverer way??
cpMAD='cp '+options.rpath+'/job.iterative.correction.madx '+options.path
system(cpMAD)

# getting madx env
#filef=os.popen("which madx")
#ddd=filef.readline()
#mad=ddd.replace("\n","")


# Line to run MAD
runMAD=options.mad+'/madx < '+options.path+'/job.iterative.correction.madx'
#runMAD='madx < ./job.iterative.correction.madx'
print runMAD

listvar=options.var.split(",")
print listvar
varslist=[]
for var in listvar:
    
    exec('variable='+var+'()')
    varslist=varslist+variable


# Iteration loop
nlimit=int(options.Iteration) # No. of iterations

for nit in range(1,nlimit+1):
    system(runMAD)
    phasesubtract(prex,prey,nit,options,options.path)
    runcorrection(nit,options,varlist)
    cpfile='cp '+options.path+'/twiss.corrected.dat '+options.path+'/twiss.'+str(nit)+'.dat'
    system(cpfile)

# Output twiss.corrected.dat at the end of iteration
system(runMAD)
mvfile='mv '+options.path+'/twiss.corrected.dat '+options.path+'/twiss.final.dat'
system(mvfile)

# Remove MAD script
rmfile='rm '+options.path+'/job.iterative.correction.madx'
system(rmfile)
#system('rm ./out.mou')

# Restore getphasex.out and getphasey.out
command='cp '+options.path+'/getphasex.00.out '+options.path+'/getphasex.out'
system(command)
command='cp '+options.path+'/getphasey.00.out '+options.path+'/getphasey.out'
system(command)

# read file and write it to a tfs table
fpknobs=options.path+'/changeparameters_all'
preknobs=open(fpknobs,'r')
predelta=[]
kq=[]
for line in preknobs:
    if line=='return;':
        break
    lines=line.split()
    kq.append(lines[0])
    predelta.append(lines[5])

preknobs.close()



tfsknobs=options.path+'/changeparameters_all.tfs'
tfsknobs=open(tfsknobs,'w')
tfsknobs.write("* NAME DELTA \n")
tfsknobs.write("$ %s %le \n")

for line in range(len(kq)):

   
    tfsknobs.write(kq[line]+" "+predelta[line]+"\n")
    
tfsknobs.close()


# change the values to LSA and YASP
if options.ACCEL=="SPS":
      
	b=twiss(options.path+"/changeparameters_all.tfs")
	execfile(options.rpath+'/Bumps.py')    # LOADS corrs
	execfile(options.rpath+'/BumpsYASP.py') # LOADS corrsYASP
        #Output for YASP...
	f=open(options.path+"/changeparameters_all.yasp", "w")
        #Output for Knob...
        g=open(options.path+"/changeparameters_all.knob", "w")
        f.write("#PLANE H\n")
	f.write("#UNIT RAD\n")
	
        g.write("* NAME  DELTA \n")
        g.write("$ %s    %le   \n")
	plane = 'H'
	beam = '1'
	for corr in corrsYASP:
		print >>f, "#SETTING", corr,  corrsYASP[corr]
	for corr in corrs:
                print >>g, "K"+corr, corrs[corr]
	f.close()
        g.close()
