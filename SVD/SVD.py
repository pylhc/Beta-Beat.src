############################################################
#
# Adjusted SVD
#
#
#
#
#
#
############################################################

# imports
from optparse import OptionParser
import sys,os
from numpy import linalg
from numpy.linalg import svd as singvd
from numpy import *
from pylab import *



# dealing with import
def parsse():
    parser = OptionParser()
    parser.add_option("-p", "--path", 
                      help="Path",
                      metavar="path", default="./",dest="path")    
    (opt, args) = parser.parse_args()

    return opt

# reading file
def readFile(inputfile):
    

    maH=[];maV=[];title=[];bpmh=[];bpmv=[];#bpmh={};bpmv={}

    counth=0;countv=0

    ifile=open(inputfile,"r")

    lines=ifile.readlines()

    for line in lines:
        char=line.split()

        if char[0]=="0": # hor
            lst=char[3:]
            lst=[float(i) for i in lst]
            maH.append(lst)
            #bpmh[char[1]]=counth
            bpmh.append(char[1])
            counth=counth+1
        elif char[0]=="1": # ver
            lst=char[3:]
            lst=[float(i) for i in lst]
            maV.append(lst)
            #bpmv[char[1]]=countv
            bpmv.append(char[1])
            countv=countv+1
        elif char[0]=="#": # titel
            title.append(line)
        else:
            print "ERROR => Unknown input",char[0]
            sys.exit()

    ifile.close()

    if len(maH)==0.0:
            print "ERROR => No turns in data file"
            sys.exit()

    return maH,maV,title,bpmh,bpmv

def readDri(path):

    drinp=path+"/Drive.inp"
    driterms=path+"/DrivingTerms"

    global kick,tunex,tuney,istun,picks,picke,eturn,inputfile
    
    if os.path.exists(drinp) and os.path.exists(driterms):

        inp=open(drinp,"r").readlines()

        for line in inp:

            char=line.split("=")

            if "KICK "==char[0]:
                kick=int(char[1])

            elif "TUNE X" in char[0]:
                tunex=float(char[1])

            elif "TUNE Y" in char[0]:
                tuney=float(char[1])
            elif "ISTUN" in char[0]:
                istun=float(char[1])

            elif "PICKUP START" in char[0]:
                picks=float(char[1])

            elif "PICKUP END" in char[0]:
                picke=float(char[1])


        term=open(driterms,"r").readlines()
        char=term[0].split()
        eturn=int(char[2])
        inputfile=char[0]
            
    else:

        print  "ERROR => Drive.inp or/and Drivingterms not found!"
        sys.exit()


def reshape(maH,maV):



   if (kick<eturn) and (eturn<=len(maH[0])):

       maH=shaper(maH,kick,eturn)
       maV=shaper(maV,kick,eturn)      

   else:

       print "WARN => ignoring kick and number turn input, check above message"
       maH=maH
       maV=maV

   return maH,maV

def shaper(ma,start,end):
    temp=[]
    for i in range(len(ma)):
        temp.append(ma[i][start:end])

    return temp


def getCO(ma):
    co=[]
    maN=[]

    for i in range(len(ma)):
        tbt=array(ma[i])
        ave=average(tbt)
        co.append(ave)
        maN.append(tbt-ave)
    
    return maN,co

def transposeMa(ma,tt):

    ma=transpose(ma)
    
    return ma/sqrt(tt)
    

def Performsvd(ma):
    
    u,s,v=singvd(ma)

    #plot(u[0])
    #plot(u[8])
    #show()
    

    return u,s,transpose(v)

def getphase(up,down):

    phase=arctan2(up,down)/2./pi
    

    return array(phase)

def getfft(U):

    d=peak(abs(fft(U[0],axis=0)))

    print d
    
    


######## main part
print "INFO => SVD, parsing information"
opt=parsse()
path=opt.path

readDri(path)
print "INFO => input file: ",inputfile," reading into memory"
maH,maV,title,bpmh,bpmv=readFile(inputfile)

print "INFO => number of turns to analyse: ",eturn,"\n Starting from turn number : ",kick,"\n Number of turns acquired",len(maH[0]),"\n Number of BPMs acquired",len(maH)
#maH,maV=reshape(maH,maV)
print "INFO => After reshaping\n number of turns to analyse: ",eturn,"\n Starting from turn number : ",kick,"\n Number of turns acquired",len(maH[0]),"\n Number of BPMs acquired",len(maH)

print "INFO => reducing to zero orbit"
maH,cox=getCO(maH)
maV,coy=getCO(maV)
print "INFO => Largest CO x:",max(abs(array(cox)))," CO y:",max(abs(array(coy)))
print "INFO => Reshaping for SVD"
maH=transposeMa(maH,eturn)
maV=transposeMa(maV,eturn)

print "INFO => Performing SVD"
Ux,Sx,Vx=Performsvd(maH)
#Uy,Sy,Vty=Performsvd(maV)

phase=getphase(Sx[0]*Vx[:,0],Sx[1]*Vx[:,1])


da=open("data.out","w")
for i in range(len(phase)):

    print >> da,bpmh[i],phase[i]

da.close()

getfft(Ux)
    

print "INFO => Getting phase"
