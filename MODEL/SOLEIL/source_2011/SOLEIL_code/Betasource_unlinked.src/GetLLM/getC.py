from metaclass import *
from optparse import OptionParser
from linreg import linreg
from math import sqrt, atan, atan2, pi
import sys

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
		result.append((x0.S[x0.indx[bpm]], bpm))
		
	result.sort()
	return result




parser = OptionParser()
#parser.add_option("-a", "--accel",
#                help="Which accelerator: LHCB1 LHCB2 SPS RHIC",
#                metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-m", "--model",
                help="Twiss File",
                metavar="TwissFile", default="0", dest="Twiss")
parser.add_option("-f", "--files",
                help="paths to the results dirs separated by comma",
                metavar="FILES", default="0", dest="files")
parser.add_option("-p", "--dpps",
                help="dp/p's for the results dirs separated by comma and in the same order",
                metavar="DPPS", default="0", dest="dpps")
parser.add_option("-u", "--dpunit",
                help="dpp unit, defaul 0.001",
                metavar="DPUNIT", default="0.001", dest="dpunit")
parser.add_option("-o", "--output",
                help="output path, defaul ./",
                metavar="OUT", default="./", dest="out")
parser.add_option("-t", "--tunes",
                help="tunes Qxd,Qx,Qyd,Qy",
                  metavar="tune", default="0,0,0,0", dest="tune")


(options, args) = parser.parse_args()


m=twiss(options.Twiss)
dpps=options.dpps.split(",")
for i in range(len(dpps)):
    dpps[i]=float(dpps[i])*float(options.dpunit)
files=options.files.split(",")


m.chrombeat()

tfilesc=[]


for el in files:
    tfilesc.append(twiss(el+"/getcouple.out"))


zeroi=-1
for i in range(len(dpps)):
    if float(dpps[i])==0.0:
        zeroi=i
        tzerox=tfilesc[i]
        #tzeroy=tfilesy[i]
        print "Found dpp=0 for case:", files[i]
if zeroi<0:
    print "dpp=0 not found, better stop"
    sys.exit()


f=open(options.out+"C.out","w")
print >>f,"* NAME","S","dcc","dccerr","dcs","dcserr"
print >>f, "$ %s  %le  %le  %le %le %le"

g=open(options.out+"Cfree.out","w")
print >>g,"* NAME","S","dcc","dccerr","dcs","dcserr", "DQmin_a", "DQminerr_a", "DQmin", "DQminerr"
print >>g, "$ %s  %le  %le  %le %le %le %le %le %le %le"

bpms=intersect(tfilesc)


#convert to get rid of ac-dipole
tunes=options.tune.split(",")
qxd=float(tunes[0])
qx=float(tunes[1])
qyd=float(tunes[2])
qy=float(tunes[3])



upac=sqrt(sin(pi*abs(qxd-qy))-sin(pi*abs(qx-qxd)))
downac=sin(pi*(qx-qy))

upq=cos(2*pi*qx)-cos(2*pi*qy)
downq=pi*(sin(2*pi*qx)+sin(2*pi*qy))




for bpm in bpms:
    el=bpm[1]
    indx=[]

    c=[]
    cs=[]

    for file in tfilesc:
        ix=file.indx[el]
        indx.append(ix)

        c.append(file.F1001W[ix])
        cs.append(file.F1010W[ix])


    #print len(dpps),len(c)
    #print dpps
    print dpps,c
    cfit=linreg(dpps,c)
    csfit=linreg(dpps,cs)
    c=cfit[0]
    cerr=cfit[3]
    cs=csfit[0]
    cserr=csfit[3]
    

    print >>f,file.NAME[ix],file.S[ix],c,cerr,cs,cserr

    
    f1001free=upac/downac*c
    f1010free=upac/downac*cs
    f1001freeE=upac/downac*cerr
    f1010freeE=upac/downac*cserr

    Qup=cos(2*pi*qx)-cos(2*pi*qy)
    Qdown=pi*(sin(2*pi*qx)+sin(2*pi*qy))

    Qmin_approx=(Qup/Qdown)*4*c
    Qminerr_approx=(Qup/Qdown)*4*c

    #print f1001free,f1010free

    try:
        rdt_fac=(4*sqrt(f1001free**2-f1010free**2))/(1+4*(f1001free**2-f1010free**2))
    except:
        rdt_fac=0

    print rdt_fac,f1001free,c

    Qmin=(Qup/Qdown)*rdt_fac
    Qminerr=(Qup/Qdown)*rdt_fac

    print >>g,file.NAME[ix],file.S[ix],f1001free,f1001freeE,f1010free,f1010freeE,Qmin_approx,Qminerr_approx,Qmin,Qminerr
    


f.close()
