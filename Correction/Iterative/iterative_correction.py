"""
IMPORTANT

iterative_correction is not used for several years. It is currently not required and thus no refactoring and improving
will be applied to this script.

It is not working with the current GUI version. It crashes on a lot of places.
If this script will be needed again one have to put effort to get it running again.
TODOs for the script would be:
 - Insert error handling to get it running
 - Clean import section and use __init__ to add root of Beta-Beat.src(see GetLLM.py for example)
 - Put global executions in a proper main() function
 - Clean the awful thing with options.JustOneBeam --> if j=='1' could then be replaced through:
     use_two_beams = options.JustOneBeam == '1'
     ...
     if use_two_beams:
         ...
 - Clean also correctTESTFUNC5.py
 - Remove versionnumbers in script files and delete old versions. Git does this job for us(pretty well)

--vimaier

"""




### Module for Iterative correction
### V1.0 created by Masa. Aiba 15/Mar/2009
### V1.1 created by Glenn Vanbavinckhove 15/June/2009
### V1.2 use correct.iterative.py to avoid the problem with the new feature in GUI, selecting the family
### V1.4 enables two beams correction
### V1.42 bug for the dispersion is fixed.
### V1.43 4/Mar/2010

import pickle
import sys
import os
import json


new_path = os.path.abspath(os.path.join( os.path.dirname(os.path.abspath(__file__)), "..", "..", "Python_Classes4MAD" ))
if new_path not in sys.path:
    sys.path.append(new_path)

from metaclass import twiss
try:
	from Numeric import *
except:
	from numpy import *
from math import *
import cmath
from os import system
#import operator
from string import *
#from AllLists import *
from optparse import OptionParser
import correct


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



def phasesubtract(prex,prey,nit,options,maccel,mpath,fmodel):

    global dictionary

    print mpath



    twissfile=mpath+"/twiss.corrected.dat"
    print "Where to read the twiss file "+twissfile
    MADTwiss1=twiss(twissfile,dictionary)
    [phase1x,phase1y]=GetPhases(prex,prey,MADTwiss1)
    MADTwiss0=twiss(fmodel,dictionary)
    print "Reading model from ",fmodel
    [phase0x,phase0y]=GetPhases(prex,prey,MADTwiss0)
    fpresentx=mpath+'/getphasex.0'+str(nit)+'.out'
    fpresenty=mpath+'/getphasey.0'+str(nit)+'.out'
    faftx=open(fpresentx,'w')
    tunex=prex.Q1-(MADTwiss1.Q1 %1)+(MADTwiss0.Q1 %1)
    if tunex>0.5: tunex=tunex-1.0
    tuney=prey.Q2-(MADTwiss1.Q2 %1)+(MADTwiss0.Q2 %1)
    if tuney>0.5: tuney=tuney-1.0
    faftx.write('@ Q1 %le '+str(tunex)+'\n')
    faftx.write('* NAME   NAME2  S      PHASEX  STDPHX  PHXMDL\n')
    faftx.write('$ %s     %s     %le    %le    %le    %le\n')
    fafty=open(fpresenty,'w')
    fafty.write('@ Q2 %le '+str(tuney)+'\n')
    fafty.write('* NAME   NAME2  S      PHASEY  STDPHY  PHYMDL\n')
    fafty.write('$ %s     %s     %le    %le    %le    %le\n')

    for i in prex.NAME:
        i2=prex.NAME2[prex.indx[i]]
        try:
            phasex=prex.PHASEX[prex.indx[i]]-(phase1x[i]-phase0x[i])
            faftx.write('"'+i+'" '+'"'+i2+'" '+str(MADTwiss0.S[MADTwiss0.indx[i]])+' '+str(phasex)+' '+str(prex.STDPHX[prex.indx[i]])+' '+str(phase0x[i])+'\n')
        except:
            0.0
    faftx.close()
    system('cp '+fpresentx+' '+mpath+'/getphasex.out.copy')

    for i in prey.NAME:
        i2=prey.NAME2[prey.indx[i]]
        #print i, MADTwiss0.S[MADTwiss0.indx[i]], prex.POS1[prex.indx[i]]
        try:
            phasey=prey.PHASEY[prey.indx[i]]-(phase1y[i]-phase0y[i])
            #print 'phase difference',(phase1y[i]-phase0y[i])
            #print (phase1[i][1]-phase0[i][1])
            fafty.write('"'+i+'" '+'"'+i2+'" '+str(MADTwiss0.S[MADTwiss0.indx[i]])+' '+str(phasey)+' '+str(prey.STDPHY[prey.indx[i]])+' '+str(phase0y[i])+'\n')
        except:
            0.0

    fafty.close()
    system('cp '+fpresenty+' '+mpath+'/getphasey.out.copy')


def dispsubtract(predx,nit,options,maccel,mpath,fmodel):

    global dictionary

    twissfile=mpath+"/twiss.corrected.dat"
    MADTwiss1=twiss(twissfile,dictionary)
    MADTwiss0=twiss(fmodel,dictionary)
    fpresentDx=mpath+'/getNDx.0'+str(nit)+'.out'
    faftdx=open(fpresentDx,'w')
    faftdx.write('* NAME   S      NDX  STDNDX  NDXMDL\n')
    faftdx.write('$ %s     %le    %le    %le    %le\n')

    for i in predx.NAME:
        #try:
        sDX=predx.NDX[predx.indx[i]]-(MADTwiss1.DX[predx.indx[i]]/sqrt(MADTwiss1.BETX[predx.indx[i]])-MADTwiss0.DX[predx.indx[i]]/sqrt(MADTwiss0.BETX[predx.indx[i]]))
        faftdx.write('"'+i+'" '+str(MADTwiss0.S[MADTwiss0.indx[i]])+' '+str(sDX)+' '+str(predx.STDNDX[predx.indx[i]])+' '+str(predx.NDXMDL[predx.indx[i]])+'\n')
        #except:
        #    0.0
    faftdx.close()
    system('cp '+fpresentDx+' '+mpath+'/getNDx.out.copy')


def runcorrection(nit,options,varlist):

    correct(options,args)

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

    #system('cp '+fpknobs+' ./')
    system('cp '+fpknobs+' '+fpknobs+'_'+str(nit))

    if options.JustOneBeam=='1':
        system('cp '+fpknobs+' '+options.path2+'/changeparameters_all')








########### START ###############


## Opions are basically the same to that of correct.py
## and are passed to correct.py.

parser = OptionParser()
parser.add_option("-a", "--accel",
		 help="What accelerator: LHCB1 LHCB2 SPS RHIC",
		 metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-b", "--accel2",
		 help="What accelerator for another beam: LHCB1 LHCB2",
		 metavar="ACCEL2", default="LHCB2",dest="ACCEL2")
parser.add_option("-t", "--tech",
		 help="Which algorithm: SVD MICADO",
		 metavar="TECH", default="SVD",dest="TECH")
parser.add_option("-n", "--ncorr",
		 help="Number of Correctors for MICADO",
		 metavar="NCORR", default=5,dest="ncorr")
parser.add_option("-p", "--path",
		 help="Path to experimental files",
		 metavar="PATH", default="./",dest="path")
parser.add_option("-q", "--path2",
		 help="Path to experimental files for the other beam",
		 metavar="PATH2", default="./",dest="path2")
parser.add_option("-c", "--cut",
                  help="Singular value cut for the generalized inverse",
                  metavar="CUT", default=0.1 , dest="cut")
parser.add_option("-e", "--errorcut",
                  help="Maximum error allowed for the phase measurement",
                  metavar="ERRORCUT", default='0.013,0.2' , dest="errorcut")
parser.add_option("-m", "--modelcut",
                  help="Maximum difference allowed between model and measured phase",
                  metavar="MODELCUT", default='0.02,0.2' , dest="modelcut")
parser.add_option("-r", "--rpath",
                  help="Path to BetaBeat repository (default is the afs repository)",
                  metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
#parser.add_option("-l", "--loc",
#                  help="location to beta beat",
#                  metavar="LOC", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="loc")
parser.add_option("-o", "--optics",
                  help="optics",
                  metavar="OPT", default="nominal.opt" , dest="OPT")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction (default is 1e-7)",
                  metavar="MinStr", default=0.0000001 , dest="MinStr")
parser.add_option("-j", "--JustOneBeam",
                  help="0 only 1 beam corrrection, 1 for 2 beams correction (default is 0)",
                  metavar="JUSTONEBEAM", default='0' , dest="JustOneBeam")
parser.add_option("-x", "--python",
                  help="Which python to run",
                  metavar="PY", default="python" , dest="pyth")
parser.add_option("-v", "--Variables",
                  help="variables split with ,",
                  metavar="var", default="MQTb1" , dest="var")

###  Additional option for the iterative method

parser.add_option("-i", "--Iteration",
                  help="Number of iteration (default is 5)",
                  metavar="ITERATION", default=5 , dest="Iteration")
parser.add_option("-g", "--BigError",
                  help="Big errors (known by SBS approach) in MAD format file, e.g., kqtl7.r3b2 = kqtl7.r3b2 + ( +0.004 ); Note that the format is fixed, i. e., the sixth culumn should be dk.",
                  metavar="BigError", default="0" , dest="BE")
parser.add_option("-d", "--dictionary",
                help="File with the BPM dictionary",
                metavar="DICT", default="0", dest="dict")
parser.add_option("-z", "--mad",
                help="MAD executable",
                metavar="MAD", default="/afs/cern.ch/group/si/slap/bin/", dest="mad")
parser.add_option("-w", "--weight",
                help="Weighting factor (phasex, phasey, betax, betay, dispersion, tunes)",
                metavar="WGT", default="1,1,0,0,0,10", dest="WGT")



(options, args) = parser.parse_args()


if options.Iteration=='0':
    print "Number of iteration cannot be equal to ZERO !! \n System will exit"
    sys.exit()


if options.dict!='0':
    execfile(options.dict)
else:
    dictionary={}


j=int(options.JustOneBeam)

#### write directory file

filename=options.path+'/var4mad.sh'
file4nad=open(filename,'w')
file4nad.write('sed    -e \'s/%filedes/\'\"'+str(options.path.replace('/','\/'))+'\"\'/g\' \\\n')
file4nad.write('<'+options.rpath+'/MODEL/LHCB/fullresponse/'+options.ACCEL+'/'+'job.iterative.correction.mask > '+options.path+'/job.iterative.correction.madx \n')
#print options.path+'/job.iterative.correction.madx '
file4nad.close()
os.system("chmod 777 "+str(filename))
os.system(str(filename))

if j==1:
    filename=options.path2+'/var4mad.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed    -e \'s/%filedes/\'\"'+str(options.path2.replace('/','\/'))+'\"\'/g\' \\\n')
    file4nad.write('<'+options.rpath+'/MODEL/LHCB/fullresponse/'+options.ACCEL2+'/'+'job.iterative.correction.mask > '+options.path2+'/job.iterative.correction.madx \n')
    file4nad.close()
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))


#sys.exit()

# Save getphasex/y.out, getDx.out for correction
if os.path.exists(options.path+'/getphasex_free.out'):
	command='cp '+options.path+'/getphasex_free.out '+options.path+'/getphasex_free.out.copy'
	command2='cp '+options.path+'/getphasex_free.out.copy '+options.path+'/getphasex.00.out'
else:
	command='cp '+options.path+'/getphasex.out '+options.path+'/getphasex.out.copy'
	command2='cp '+options.path+'/getphasex.out.copy '+options.path+'/getphasex.00.out'

try:
    system(command)
    system(command2)
except:
    print 'No getphasex.out file in ',options.path
    print 'Correction is not possible... leaving Iterative Correction.'
    sys.exit()

if j==1:
    if os.path.exists(options.path2+'/getphasex_free.out'):
	    command='cp '+options.path2+'/getphasex_free.out '+options.path2+'/getphasex.out.copy'
	    command2='cp '+options.path2+'/getphasex.out.copy '+options.path2+'/getphasex.00.out'
    else:
	    command='cp '+options.path2+'/getphasex.out '+options.path2+'/getphasex.out.copy'
	    command2='cp '+options.path2+'/getphasex.out.copy '+options.path2+'/getphasex.00.out'
    try:
        system(command)
        system(command2)
    except:
        print 'No getphasex.out file in ',options.path2
        print 'Correction is not possible... leaving Iterative Correction.'
        sys.exit()

if os.path.exists(options.path2+'/getphasex_free.out'):
	command='cp '+options.path+'/getphasey_free.out '+options.path+'/getphasey.out.copy'
	command2='cp '+options.path+'/getphasey.out.copy '+options.path+'/getphasey.00.out'
else:
	command='cp '+options.path+'/getphasey.out '+options.path+'/getphasey.out.copy'
	command2='cp '+options.path+'/getphasey.out.copy '+options.path+'/getphasey.00.out'

try:
    system(command)
    system(command2)
except:
    print 'No getphasey.out file in ',options.path
    print 'Correction is not possible... leaving Iterative Correction.'
    sys.exit()

if j==1:
    if os.path.exists(options.path2+'/getphasex_free.out'):
	command='cp '+options.path2+'/getphasey_free.out '+options.path2+'/getphasey.out.copy'
	command2='cp '+options.path2+'/getphasey.out.copy '+options.path2+'/getphasey.00.out'
    else:
	command='cp '+options.path2+'/getphasey.out '+options.path2+'/getphasey.out.copy'
	command2='cp '+options.path2+'/getphasey.out.copy '+options.path2+'/getphasey.00.out'
    try:
        system(command)
        system(command2)
    except:
        print 'No getphasey.out file in ',options.path2
        print 'Correction is not possible... leaving Iterative Correction.'
        sys.exit()



#wi=int(options.WGT.split(",")[4]) # Fractional must be accepted but 1 or 0 at this moment.
wi=float(options.WGT.split(",")[4]) # Fractional must be accepted but 1 or 0 at this moment.
command='cp '+options.path+'/getNDx.out '+options.path+'/getNDx.out.copy'
command2='cp '+options.path+'/getNDx.out.copy '+options.path+'/getNDx.00.out'
try:
    fndx=options.path+'/getNDx.out'
    ftemp=open(fndx,'r')
    ftemp.close() # redundant check
    system(command)
    system(command2)
    fragndx=1
except:
    print 'No getNDx.out file in ',options.path
    print 'We continue correction ignoring the dispersion.'
    fragndx=0
if fragndx==0 and wi==1:
    print 'You select the option (or in default) to correct the dispersion, but there seems no measurement.'
if fragndx==1 and wi==0:
    print 'You have the dispersion measurement, but you are not going to correct.'
    fragndx=0

if j==1:
    command='cp '+options.path2+'/getNDx.out '+options.path2+'/getNDx.out.copy'
    command2='cp '+options.path2+'/getNDx.out.copy '+options.path2+'/getNDx.00.out'
    try:
        fndx=options.path2+'/getNDx.out'
        ftemp=open(fndx,'r')
        ftemp.close()
        system(command)
        system(command2)
        fragndx=1
    except:
        print 'No getNDx.out file in ',options.path2
        print 'We continue correction ignoring the dispersion.'
        fragndx=0
    if fragndx==0 and wi==1:
        print 'You select the option (or in default) to correct the dispersion, but there seems no measurement.'
    if fragndx==1 and wi==0:
        print 'You have the dispersion measurement, but you are not going to correct.'
        fragndx=0





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
system('cp '+file1+' '+options.path+'/changeparameters_all_0')
if j==1:
    system('cp '+file1+' '+options.path2+'/changeparameters_all')



# Load measurement data
filex=options.path+'/getphasex.00.out'
prex=twiss(filex)
filey=options.path+'/getphasey.00.out'
prey=twiss(filey)
if fragndx==1:
    try:
        filedx=options.path+'/getNDx.00.out'
        predx=twiss(filedx)
    except:
        print "getNDx not found"
        fragndx=0

if j==1:
    filex2=options.path2+'/getphasex.00.out'
    prex2=twiss(filex2)
    filey2=options.path2+'/getphasey.00.out'
    prey2=twiss(filey2)
    if fragndx==1:
        try:
            filedx2=options.path2+'/getNDx.00.out'
            predx2=twiss(filedx2)
        except:
            print "getNDx not found"
            fragndx=0




# Copy MAD script to the current directory temporarily since MAD output are generated in the currect directory.
# This avoids to output files in the repository.
# Is there any other cleverer way??
#cpMAD='cp '+options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/job.iterative.correction.madx '+options.path
#system(cpMAD)

# getting madx env
#filef=os.popen("which madx")
#ddd=filef.readline()
#mad=ddd.replace("\n","")


# Line to run MAD
runMAD=options.mad+' < '+options.path+'/job.iterative.correction.madx > out.mou'
if j==1:
    runMAD2=options.mad+' < '+options.path2+'/job.iterative.correction.madx > out.mou'
# import knobs list
path_to_json_file = os.path.join(options.rpath, "MODEL", "LHCB", "fullresponse", options.ACCEL, "AllLists.json")
knobsdict=json.load(file(path_to_json_file,'r'))
if j==1:
    path_to_json_file = os.path.join(options.rpath, "MODEL", "LHCB", "fullresponse", options.ACCEL2, "AllLists.json")
    knobsdict=json.load(file(path_to_json_file,'r'))
    #TODO: overrides knobsdict reference of first options.ACCEL. Seems not to be correct.
    # AllLists have same keys but different values. A union of both dicts seems to be the solution
    # -- vimaier


listvar=options.var.split(",")
print listvar
varlist=[]
for var in listvar:
    variable=knobsdict[var]
    varlist=varlist+variable


# Iteration loop
nlimit=int(options.Iteration) # No. of iterations

for nit in range(1,nlimit+1):
    print "Sending command to madx :",runMAD
    system(runMAD)
    bn1=upper(prex.NAME[0])
    bn2=upper(prex.NAME2[0])
    print bn1,bn2
    fmodel=options.OPT.split(",")[0]+"/twiss.dat"
    print fmodel
    phasesubtract(prex,prey,nit,options,options.ACCEL,options.path,fmodel)
    if fragndx==1:
        dispsubtract(predx,nit,options,options.ACCEL,options.path,fmodel)
    if j==1:
        system(runMAD2)
	print "Running mad2"
	fmodel=options.OPT.split(",")[1]+"/twiss.dat"
        phasesubtract(prex2,prey2,nit,options,options.ACCEL2,options.path2,fmodel)
        if fragndx==1:
            dispsubtract(predx2,nit,options,options.ACCEL2,options.path2,fmodel)
    runcorrection(nit,options,varlist)
    cpfile='cp '+options.path+'/twiss.corrected.dat '+options.path+'/twiss.'+str(nit)+'.dat'
    system(cpfile)
    if j==1:
        cpfile='cp '+options.path2+'/twiss.corrected.dat '+options.path2+'/twiss.'+str(nit)+'.dat'
        system(cpfile)


# Output twiss.corrected.dat at the end of iteration
system(runMAD)
mvfile='mv '+options.path+'/twiss.corrected.dat '+options.path+'/twiss.final.dat'
system(mvfile)
if j==1:
    system(runMAD2)
    mvfile='mv '+options.path2+'/twiss.corrected.dat '+options.path2+'/twiss.final.dat'
    system(mvfile)

# Remove MAD script
#rmfile='rm '+options.path+'/job.iterative.correction.madx'
#system(rmfile)
#system('rm ./out.mou')

# Restore getphasex.out and getphasey.out, no need to restore anymore as correct.py (TESTFUNC) works on the copied files.
#command='cp '+options.path+'/getphasex.00.out '+options.path+'/getphasex.out'
#system(command)
#command='cp '+options.path+'/getphasey.00.out '+options.path+'/getphasey.out'
#system(command)
#if fragndx==1:
#    command='cp '+options.path+'/getNDx.00.out '+options.path+'/getNDx.out'
#    system(command)
#if j==1:
#    command='cp '+options.path2+'/getphasex.00.out '+options.path2+'/getphasex.out'
#    system(command)
#    command='cp '+options.path2+'/getphasey.00.out '+options.path2+'/getphasey.out'
#    system(command)
#    if fragndx==1:
#        command='cp '+options.path2+'/getNDx.00.out '+options.path2+'/getNDx.out'
#        system(command)


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
    predelta.append(-float(lines[5]))

preknobs.close()



tfsknobs=options.path+'/changeparameters_all.tfs'
tfsknobs=open(tfsknobs,'w')
tfsknobs.write("* NAME DELTA \n")
tfsknobs.write("$ %s %le \n")

for line in range(len(kq)):


    tfsknobs.write(kq[line]+" "+str(predelta[line])+"\n")

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
