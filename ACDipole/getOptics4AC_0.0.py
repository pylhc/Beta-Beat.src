################################################################
#                                                              #
#  @ Glenn Vanbavinckhove  (gvanbavi@cern.ch)=> Date: 11/02/10 #
#  @ Ryoichi Miyamoto (miyamoto@bnl.gov)                       #
#                                                              #
################################################################
#
#  !=> getOptics4AC_0.0.py : - Construction of the program (11/02/10)

# imports
from metaclass import twiss
from optparse import OptionParser
from math import *
import sys,os

# option parser
parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                metavar="ACCEL", default="LHCB1",dest="accel")
parser.add_option("-p", "--path",
                help="Path to output files of GetLLM",
                metavar="PATH", default="./",dest="path")
parser.add_option("-b", "--bbsrouce",
                help="beta beat source",
                metavar="bb", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/", dest="bb")
parser.add_option("-o", "--output path",
                help="Where to put the output files",
                metavar="output", default="./", dest="output")
parser.add_option("-d", "--drive tune",
                help="input of drive tune for horizontal and vertical",
                metavar="dtune", default="0.28,0.31", dest="dtune")

(options, args) = parser.parse_args()

#function
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


def getbeta(lambdad,phasebpm,phaseac,phasedriven,driventune,drivenbeta,sign):


    if sign=="+":
        forcos=2*(phasebpm-phaseac-phasedriven+pi*driventune)
    else:
        forcos=2*(phasebpm-phaseac-phasedriven-pi*driventune)

    up=1+lambdad**2+2*lambdad*cos(forcos)
    down=1-lambdad**2
    
    beta=(up/down)*drivenbeta


    return beta

def getalfa(lambdad,phasebpm,phaseac,phasedriven,driventune,drivenalfa,sign):

    if sign=="+":
        forcos=2*(phasebpm-phaseac-phasedriven+pi*driventune)
    else:
        forcos=2*(phasebpm-phaseac-phasedriven-pi*driventune)

    up1=1+lambdad**2+2*lambdad*cos(forcos)
    up2=2*lambdad*sin(forcos)
    down=1-lambdad**2

    alfa=(up1/down)*drivenalfa+(up2/down)

    return alfa


def getphase(sign,lambdad,phasebpm,phaseac,phasedriven,driventune):

    if sign=="+":
        if phasebpm<pi:
            up=(1-lambdad**2)*sin(phasebpm)
            down=(1+lambdad**2)*cos(phasebpm)+2*lambdad*cos(phasebpm-2*phaseac-2*phasedriven+2*pi*driventune)
            phase=atan(up/down)
        else:
            up=(1-lambdad**2)*sin(phasebpm)
            down=(1+lambdad**2)*cos(phasebpm)+2*lambdad*cos(phasebpm-2*phaseac-2*phasedriven+2*pi*driventune)
            phase=atan(up/down)+pi*(phasebpm-pi)
    else:
        if phasebpm<pi:
            up=(1-lambdad**2)*sin(phasebpm-2*pi*driventune)
            down=(1+lambdad**2)*cos(phasebpm-2*pi*driventune)+2*lambdad*cos(phasebpm-2*phaseac-2*phasedriven)
            phase=atan(up/down)
        else:
            up=(1-lambdad**2)*sin(phasebpm-2*pi*driventune)
            down=(1+lambdad**2)*cos(phasebpm-2*pi*driventune)+2*lambdad*cos(phasebpm-2*phaseac-2*phasedriven)
            phase=atan(up/down)+pi*(phasebpm-2*pi*driventune-pi)+2*pi*driventune

    return phase


def getdrivenphase(phasemodel,tune,lambdad,driventune):

    phasetune=phasemodel-pi*tune
    lambdafun=(1+lambdad)/(1-lambdad)
    
    if phasetune<pi:
        phasedriven=atan(lambdafun*tan(phasetune))
    else:
        phasedriven=atan(lambdafun*tan(phasetune))+pi*(phasetune-pi)+pi*driventune

    return phasedriven
    


def getlambda(drivetune,tune):

    up=sin(pi*(drivetune-tune))
    down=sin(pi*(drivetune+tune))
    lambdad=up/down

    return lambdad


# reading files
phasex=twiss(options.path+'/getphasex.out')
Qx=phasex.Q1
Qy=phasex.Q2
phasey=twiss(options.path+'/getphasey.out')
phasext=twiss(options.path+'/getphasetotx.out')
phaseyt=twiss(options.path+'/getphasetoty.out')
betax=twiss(options.path+'/getbetax.out')
betay=twiss(options.path+'/getbetay.out')
ampbetax=twiss(options.path+'/getampbetax.out')
ampbetay=twiss(options.path+'/getampbetay.out')

model=twiss(options.bb+'/MODEL/'+options.accel+'/nominal.opt/twiss_ac.dat')
allfiles=[phasex,phasey,phasext,phaseyt,betax,betay,ampbetax,ampbetay]

# main part
bpms=intersect(allfiles)
bpms=modelIntersect(bpms, model)


for bpm in bpms:

    # horizontal
    phasemodel=model.MUX[model.indx["HAC"]]-model.MUX[model.indx["HAC"]-1]
    drivetune=float(options.dtune.split(',')[0])
    lambdad=getlambda(drivetune,Qx)
    drivenphase=getdrivenphase(phasemodel,Qx,lambdad,drivetune)

    phasebpm=phasext.PHASEX[phasext.indx[bpm[1]]]
    phaseac=phasext.PHASEX[phasext.indx["BPMYB.6L4.B1"]]
    drivebeta=betax.BETX[betax.indx[bpm[1]]]
    drivealfa=betax.ALFX[betax.indx[bpm[1]]]

    sac=model.indx["HAC"]
    s=model.indx[bpm[1]]
    if(s<sac):
        sign="+"
    else:
        sign="-"
        
    beta=getbeta(lambdad,phasebpm,phaseac,drivenphase,drivetune,drivebeta,sign)
    alfa=getalfa(lambdad,phasebpm,phaseac,drivenphase,drivetune,drivealfa,sign)
    phase=getphase(sign,lambdad,phasebpm,phaseac,drivenphase,drivetune)


    horizontal=[bpm, beta,alfa,phase]

    # vertical
    phasemodel=model.MUY[model.indx["VAC"]]-model.MUY[model.indx["VAC"]-1]
    drivetune=float(options.dtune.split(',')[0])
    lambdad=getlambda(drivetune,Qy)
    drivenphase=getdrivenphase(phasemodel,Qy,lambdad,drivetune)

    phasebpm=phaseyt.PHASEY[phaseyt.indx[bpm[1]]]
    phaseac=phaseyt.PHASEY[phaseyt.indx["BPMYB.6L4.B1"]]
    drivebeta=betay.BETY[betay.indx[bpm[1]]]
    drivealfa=betay.ALFY[betay.indx[bpm[1]]]

    sac=model.indx["VAC"]
    s=model.indx[bpm[1]]
    if(s<sac):
        sign="+"
    else:
        sign="-"
        
    beta=getbeta(lambdad,phasebpm,phaseac,drivenphase,drivetune,drivebeta,sign)
    alfa=getalfa(lambdad,phasebpm,phaseac,drivenphase,drivetune,drivealfa,sign)
    phase=getphase(sign,lambdad,phasebpm,phaseac,drivenphase,drivetune)   

    vertical=[bpm, beta,alfa,phase]

    print vertical
