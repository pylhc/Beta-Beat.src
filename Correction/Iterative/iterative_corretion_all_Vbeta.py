### Module for Iterative correction
### V1.0 created by Masa. Aiba 15/Mar/2009

from metaclass import *
from Numeric import *
from math import *
import cmath
import sys, pickle,os
from os import system
#import operator
from string import *
#from AllLists import *
from optparse import OptionParser



#------------


def intersect(ListOfFile): 
	'''Pure intersection of all bpm names in all files '''
	if len(ListOfFile)==0:
		print "Nothing to intersect!!!!"
		sys.exit()
	z=ListOfFile[0].NAME
	for b in ListOfFile:
		z=filter(lambda x: x in z   , b.NAME)
	#SORT by S
	result=[]
	x0=ListOfFile[0]
	for bpm in z:
		try:
			result.append((x0.S[x0.indx[bpm]], bpm))
		except:
			result.append((x0.POS[x0.indx[bpm]], bpm))
		
	result.sort()
	return result

def modelIntersect(expbpms, model):
	bpmsin=[]
	for bpm in expbpms:
		try:
			check=model.indx[bpm[1].upper()]
			bpmsin.append(bpm)
		except:
			print bpm, "Not in Model"
	if len(bpmsin)==0:
		print "Zero intersection of Exp and Model"
		print "Please, provide a good Dictionary"
		print "Now we better leave!"
		sys.exit()			
	return bpmsin




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

    global dictionary


    MADTwiss1=twiss("./twiss.corrected.dat",dictionary)
    [phase1x,phase1y]=GetPhases(prex,prey,MADTwiss1)
    fmodel=options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/twiss.dat'
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
	    print 'phase difference',(phase1x[i]-phase0x[i]) 
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


##########################################################################

def Dxsubtract(pre,nit,options):

    global dictionary


    MADTwiss1=twiss("./twiss.corrected.dat",dictionary)
    fmodel=options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/twiss.dat'
    MADTwiss0=twiss(fmodel,dictionary)
    fpresent=options.path+'/getNDx.0'+str(nit)+'.out'
    faft=open(fpresent,'w')
    faft.write('* NAME   S      NDX    STDNDX  DX\n')
    faft.write('$ %s     %le    %le    %le     %le\n')

    expbpms=intersect([pre])
    bpms=modelIntersect(expbpms, MADTwiss0)

    for i in range(0,len(bpms)):
        bn1=upper(bpms[i][1])
        try:
		ndx0=MADTwiss0.DX[MADTwiss0.indx[bn1]]/sqrt(MADTwiss0.BETX[MADTwiss0.indx[bn1]])
		ndx1=MADTwiss1.DX[MADTwiss1.indx[bn1]]/sqrt(MADTwiss1.BETX[MADTwiss1.indx[bn1]])
		ndx=pre.NDX[pre.indx[bn1]]-(ndx1-ndx0)
		dx=pre.DX[pre.indx[bn1]]-(MADTwiss1.DX[MADTwiss1.indx[bn1]]-MADTwiss0.DX[MADTwiss0.indx[bn1]])
		faft.write('"'+bn1+'" '+str(MADTwiss0.S[MADTwiss0.indx[bn1]])+' '+str(ndx)+' '+str(pre.STDNDX[pre.indx[bn1]])+' '+str(dx)+'\n')
        except:
		0.0
    faft.close()
    system('cp '+fpresent+' '+options.path+'/getNDx.out')

##########################################################################

def Dysubtract(pre,nit,options):

    global dictionary


    MADTwiss1=twiss("./twiss.corrected.dat",dictionary)
    fmodel=options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/twiss.dat'
    MADTwiss0=twiss(fmodel,dictionary)
    fpresent=options.path+'/getDy.0'+str(nit)+'.out'
    faft=open(fpresent,'w')
    faft.write('* NAME   S      DY     STDDY\n')
    faft.write('$ %s     %le    %le    %le\n')

    expbpms=intersect([pre])
    bpms=modelIntersect(expbpms, MADTwiss0)

    for i in range(0,len(bpms)):
        bn1=upper(bpms[i][1])
        try:
		dy=pre.DY[pre.indx[bn1]]-(MADTwiss1.DY[MADTwiss1.indx[bn1]]-MADTwiss0.DY[MADTwiss0.indx[bn1]])
		faft.write('"'+bn1+'" '+str(MADTwiss0.S[MADTwiss0.indx[bn1]])+' '+str(dy)+' '+str(pre.STDDY[pre.indx[bn1]])+'\n')
        except:
		0.0
    faft.close()
    system('cp '+fpresent+' '+options.path+'/getDy.out')


################################################################################

def couplesubtract(pre,nit,options):


    global dictionary

    MADTwiss1=twiss("./twiss.corrected.dat",dictionary)
    fmodel=options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/twiss.dat'
    MADTwiss0=twiss(fmodel,dictionary)
    fpresent=options.path+'/getcouple.0'+str(nit)+'.out'
    faft=open(fpresent,'w')
    faft.write('* NAME  S      F1001R F1001I F1010R F1010I\n')
    faft.write('$ %s    %le    %le    %le    %le    %le\n')

    expbpms=intersect([pre])
    bpms=modelIntersect(expbpms, MADTwiss0)
    print bpms

    MADTwiss0.Cmatrix()
    MADTwiss1.Cmatrix()

    for i in range(0,len(bpms)):
        bn1=upper(bpms[i][1])
        try:
		f1001r=pre.F1001R[pre.indx[bn1]]-(MADTwiss1.F1001R[MADTwiss1.indx[bn1]]-MADTwiss0.F1001R[MADTwiss0.indx[bn1]])
		f1001i=pre.F1001I[pre.indx[bn1]]-(MADTwiss1.F1001I[MADTwiss1.indx[bn1]]-MADTwiss0.F1001I[MADTwiss0.indx[bn1]])
		f1010r=pre.F1010R[pre.indx[bn1]]-(MADTwiss1.F1010R[MADTwiss1.indx[bn1]]-MADTwiss0.F1010R[MADTwiss0.indx[bn1]])
		f1010i=pre.F1010I[pre.indx[bn1]]-(MADTwiss1.F1010I[MADTwiss1.indx[bn1]]-MADTwiss0.F1010I[MADTwiss0.indx[bn1]])
		faft.write('"'+bn1+'" '+str(MADTwiss0.S[MADTwiss0.indx[bn1]])+' '+' '+str(f1001r)+' '+str(f1001i)+' '+' '+str(f1010r)+' '+str(f1010i)+' '+'\n')
        except:
		0.0
    faft.close()
    system('cp '+fpresent+' '+options.path+'/getcouple.out')


#########################################################################################################


def runcorrection(nit,options,varall,fbfact):

    # Not clever implementation obliged by the fact that execfile("script.py") seems not to accept options for script.py.
    optlist=[['-a ',options.ACCEL],['-o ',options.OPT],['-t ',options.TECH],['-n ',options.ncorr],['-p ',options.path],['-c ',options.cut],['-e ',options.errorcut],['-m ',options.modelcut],['-r ',options.rpath],['-s ',options.MinStr],['-j ',options.JustOneBeam]]
    #runcorrect='pythonafs '+options.rpath+'/Correction/correct.py'
    runcorrect='pythonafs ./correct.py'
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
        delta[str(kq)]=str(float(lines[5])*fbfact) # not a 100% correction
    knobs.close()


    fpknobs=options.path+'/changeparameters_all'
    rmpknobs='rm '+fpknobs
    system(rmpknobs)

    
    newknobs=open(fpknobs,'w')
    for i in range(0,len(varall)):
        try:
            predel=float(predelta[varall[i]])
        except:
            predel=0.0
        try:
            curdel=float(delta[varall[i]])
        except:
            curdel=0.0
        delall=predel-curdel # inverse correction
        newknobs.write(varall[i]+' = '+varall[i]+' + ( '+str(delall)+' );\n')
    newknobs.write('return;')
    newknobs.close()

    system('cp '+fpknobs+' ./')
    system('cp '+file1+' '+file1+'_'+str(nit))


#########################################################################################################


def runcorrectcouple(nit,options,varall,fbfact):

    # Not clever implementation obliged by the fact that execfile("script.py") seems not to accept options for script.py.
    optlist=[['-a ',options.ACCEL],['-o ',options.OPT],['-d ',options.Dy],['-p ',options.path],['-c ',options.cut],['-r ',options.rpath]]
    #runcorrect='pythonafs '+options.rpath+'/Correction/correct_coupleDy.py'
    runcorrect='pythonafs ./correct_coupleDy.py'
    for i in range(0,len(optlist)):
        runcorrect=runcorrect+' '+optlist[i][0]+str(optlist[i][1])
    system(runcorrect)
    print runcorrect

    system('cp changeparameters changeparameters.save')
    system('cp changeparameters_couple changeparameters')

    
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
        delta[str(kq)]=str(float(lines[5])*fbfact) # not a 100% correction
    knobs.close()


    fpknobs=options.path+'/changeparameters_all'
    rmpknobs='rm '+fpknobs
    system(rmpknobs)


    newknobs=open(fpknobs,'w')
    for i in range(0,len(varall)):
        try:
            predel=float(predelta[varall[i]])
        except:
            predel=0.0
        try:
            curdel=float(delta[varall[i]])
        except:
            curdel=0.0
        delall=predel-curdel # inverse correction
        newknobs.write(varall[i]+' = '+varall[i]+' + ( '+str(delall)+' );\n')
    newknobs.write('return;')
    newknobs.close()

    system('cp '+fpknobs+' ./')
    system('cp '+file1+' '+file1+'_'+str(nit))





########### START ###############


## Opions are basically the same to that of correct.py
## and are passed to correct.py (and for correct_couple.py)

parser = OptionParser()
parser.add_option("-a", "--accel", 
		 help="What accelerator: LHCB1 LHCB2 SPS RHIC",
		 metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-o", "--opt", 
		 help="Which optics: Injection, Collision etc",
		 metavar="OPT", default="/",dest="OPT")
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
                  metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction (default is 1e-6)",
                  metavar="MinStr", default=0.000001 , dest="MinStr")
parser.add_option("-j", "--JustOneBeam",
                  help="0 Just quads from one beam are used, 1 all quads (default is 0)",
                  metavar="JUSTONEBEAM", default=0 , dest="JustOneBeam")

## for correct_coupleDy.py
parser.add_option("-y", "--Dy",
                  help="Turn-on=1 or Turn-off=0 vertical dispersion correction",
                  metavar="Dy", default=0, dest="Dy")


###  Additional option for the iterative method

parser.add_option("-i", "--Iteration",
                  help="Number of iteration (default is 9 with the feedback factor of 0.3)",
                  metavar="ITERATION", default=9 , dest="Iteration")
parser.add_option("-b", "--BigError",
                  help="Big errors (known by SBS approach) in MAD format file, e.g., kqtl7.r3b2 = kqtl7.r3b2 + ( +0.004 ); Note that the format is fixed, i. e., the sixth culumn should be dk.",
                  metavar="BigError", default="0" , dest="BE")
parser.add_option("-d", "--dictionary",
                help="File with the BPM dictionary",
                metavar="DICT", default="0", dest="dict")


(options, args) = parser.parse_args()

## reduce correction knobs to 30% 
fbfact=0.3

if options.dict!='0':
    execfile(options.dict)
else:
    dictionary={}



# Save getphasex/y.out, getNDx.out for correction

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
    system(command)
    fragndx=1
except:
    print 'No getNDx.out file in ',options.path
    print 'We continue correction ignoring the dispersion.'




# Save getcouple.out and getDy.out

command='cp '+options.path+'/getcouple.out '+options.path+'/getcouple.00.out'
try:
    system(command)
except:
    print 'No getcouple.out file in ',options.path
    print 'Correction is not possible... leaving Iterative Correction.'
    sys.exit()

command='cp '+options.path+'/getDy.out '+options.path+'/getDy.00.out'
fragdy=0
try:
    system(command)
    fragdy=1
except:
    print 'No getDy.out file in ',options.path
    print 'Correction is not possible... leaving Iterative Correction.'
    sys.exit()



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

# Initialize coupling knobs
#file1=options.path+'/changeparameters_couple_all'
#iniknobs=open(file1,'w')
#iniknobs.write('return;')
#iniknobs.close()
#system('cp '+file1+' ./')
#system('cp '+file1+' ./changeparameters_0')
#system('cp '+file1+' '+file1+'_0')


# Load measurement data
filex=options.path+'/getphasex.00.out'
prex=twiss(filex)
filey=options.path+'/getphasey.00.out'
prey=twiss(filey)
if (fragndx==1):
    filendx=options.path+'/getNDx.00.out'
    prendx=twiss(filendx)
filecouple=options.path+'/getcouple.00.out'
precouple=twiss(filecouple)
if (fragdy==1):
    filedy=options.path+'/getDy.00.out'
    predy=twiss(filedy)


# Copy MAD script to the current directory temporarily since MAD output are generated in the currect directory.
# This avoids to output files in the repository.
# Is there any other cleverer way??
cpMAD='cp '+options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/job.iterative.correction.t.madx ./'
system(cpMAD)

# Line to run MAD
runMAD='madx < ./job.iterative.correction.t.madx > out.mou'
#runMAD='madx < ./job.iterative.correction.madx'

# import knobs list
execfile(options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/AllLists.py')
varname='vars'+options.ACCEL
exec('varlist='+varname+'()')

# import coupling knobs list
execfile(options.rpath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/Coupling/AllLists_couple.py')
varname='squadvars'+options.ACCEL
exec('varcouple='+varname+'()')

varall=varlist+varcouple
print varall

# Iteration loop
nlimit=int(options.Iteration) # No. of iterations

for nit in range(1,nlimit+1):
    system(runMAD)
    phasesubtract(prex,prey,nit,options)
    print 'kiteru1'
    couplesubtract(precouple,nit,options)
    print 'kiteru2'
    if (fragndx==1): Dxsubtract(prendx,nit,options)
    print 'kiteru3'
    if (fragdy==1): Dysubtract(predy,nit,options)
    print 'kiteru4'
    runcorrection(nit,options,varall,fbfact)
    runcorrectcouple(nit,options,varall,fbfact)
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
try:
	command='cp '+options.path+'/getNDx.00.out '+options.path+'/getNDx.out'
	system(command)
except:
	0.0
try:
	command='cp '+options.path+'/getDy.00.out '+options.path+'/getDy.out'
	system(command)
except:
	0.0
command='cp '+options.path+'/getcouple.00.out '+options.path+'/getcouple.out'
system(command)
