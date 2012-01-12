#!/usr/bin/env python2.5
from os import system
import sys, StringIO
from os.path import exists
from string import split,replace
import pickle, random
stone="injectionStone"
#stone="storeStone"

qTrimNames=['bi9-tq4','bo10-tq4','bi12-tq4','bo11-tq4','bi1-tq4',\
            'bo2-tq4','bi4-tq4','bo3-tq4','bi5-tq4','bo6-tq4',\
            'bi8-tq4','bo7-tq4','bi9-tq5','bo10-tq5','bi12-tq5',\
            'bo11-tq5','bi1-tq5','bo2-tq5','bi4-tq5','bo3-tq5',\
            'bi5-tq5','bo6-tq5','bi8-tq5','bo7-tq5','bo10-tq6',\
            'bi9-tq6','bi12-tq6','bo11-tq6','bi1-tq6','bo2-tq6',\
            'bi4-tq6','bo3-tq6','bi5-tq6','bo6-tq6','bi8-tq6','bo7-tq6']

quadNames='bi9-qf1'+' '+'bo10-qd1'+' '+'bi12-qf1'+' '+'bo11-qd1'+' '\
           +'bi1-qf1'+' '+'bo2-qd1'+' '+'bi4-qf1'+' '+'bo3-qd1'+' '\
           +'bi5-qf1'+' '+'bo6-qd1'+' '+'bi8-qf1'+' '+'bo7-qd1'+' '\
           +'bi9-qd2'+' '+'bo10-qf2'+' '+'bi12-qd2'+' '+'bo11-qf2'+' '\
           +'bi1-qd2'+' '+'bo2-qf2'+' '+'bi4-qd2'+' '+'bo3-qf2'+' '\
           +'bi5-qd2'+' '+'bo6-qf2'+' '+'bi8-qd2'+' '+'bo7-qf2'+' '\
           +'bi9-qf3'+' '+'bo10-qd3'+' '+'bi12-qf3'+' '+'bo11-qd3'+' '\
           +'bi1-qf3'+' '+'bo2-qd3'+' '+'bi4-qf3'+' '+'bo3-qd3'+' '\
           +'bi5-qf3'+' '+'bo6-qd3'+' '+'bi8-qf3'+' '+'bo7-qd3'+' '\
           +'bi9-qd4'+' '+'bi12-qd4'+' '+'bi1-qd4'+' '+'bi4-qd4'+' '\
           +'bi5-qd4'+' '+'bi8-qd4'+' '+'bo10-qd5'+' '+'bo3-qd5'+' '\
           +'bi9-qf7'+' '+'bi12-qf7'+' '+'bi1-qf7'+' '+'bi4-qf7'+' '\
           +'bi5-qf7'+' '+'bi8-qf7'+' '+'bo3-qd7'+' '+'bo10-qd7'+' '\
           +'bi1-qd10'+' '+'bi1-qf11'+' '+'bo10-qd9'+' '+'bi9-qd8'+' '\
           +'bi12-qd8'+' '+'bi1-qd8'+' '+'bi4-qd8'+' '+'bi5-qd8'+' '\
           +'bi8-qd8'+' '+'bi9-qf9'+' '+'bi12-qf9'+' '+'bi1-qf9'+' '\
           +'bi4-qf9'+' '+'bi5-qf9'+' '+'bi8-qf9'+' '+'bo10-qf8'+' '\
           +'bo11-qf8'+' '+'bo2-qf8'+' '+'bo3-qf8'+' '+'bo6-qf8'+' '\
           +'bo7-qf8'+' '+'bi9-tq4'+' '+'bo10-tq4'+' '+'bi12-tq4'+' '\
           +'bo11-tq4'+' '+'bi1-tq4'+' '+'bo2-tq4'+' '+'bi4-tq4'+' '\
           +'bo3-tq4'+' '+'bi5-tq4'+' '+'bo6-tq4'+' '+'bi8-tq4'+' '\
           +'bo7-tq4'+' '+'bi9-tq5'+' '+'bo10-tq5'+' '+'bi12-tq5'+' '\
           +'bo11-tq5'+' '+'bi1-tq5'+' '+'bo2-tq5'+' '+'bi4-tq5'+' '\
           +'bo3-tq5'+' '+'bi5-tq5'+' '+'bo6-tq5'+' '+'bi8-tq5'+' '\
           +'bo7-tq5'+' '+'bo10-tq6'+' '+'bi9-tq6'+' '+'bi12-tq6'+' '\
           +'bo11-tq6'+' '+'bi1-tq6'+' '+'bo2-tq6'+' '+'bi4-tq6'+' '\
           +'bo3-tq6'+' '+'bi5-tq6'+' '+'bo6-tq6'+' '+'bi8-tq6'+' '+'bo7-tq6'
bpmNamesX='G6-BX'+' '+'BO6-B1'+' '+'BO6-B3'+' '+'BO6-B3.1'+' '+'BO6-B4'\
          +' '+'BO6-BH6'+' '+'BO6-B7'+' '+'BO6-B8'+' '+'BO6-BH10'\
          +' '+'BO6-BH12'+' '+'BO6-BH14'+' '+'BO6-BH16'+' '+'BO6-BH18'\
          +' '+'BO6-BH20'+' '+'BO7-BH20'+' '+'BO7-BH18'+' '+'BO7-BH16'\
          +' '+'BO7-BH14'+' '+'BO7-BH12'+' '+'BO7-BH10'+' '+'BO7-B8'\
          +' '+'BO7-B7'+' '+'BO7-BH6'+' '+'BO7-B4'+' '+'BO7-B3.1'\
          +' '+'BO7-B3'+' '+'BO7-B1'+' '+'G7-BX'+' '+'G8-BX'+' '+'BI8-B1'\
          +' '+'BI8-B3'+' '+'BI8-B3.1'+' '+'BI8-B4'+' '+'BI8-BH5'\
          +' '+'BI8-B7'+' '+'BI8-B8'+' '+'BI8-BH9'+' '+'BI8-BH11'\
          +' '+'BI8-BH13'+' '+'BI8-BH15'+' '+'BI8-BH17'+' '+'BI8-BH19'\
          +' '+'BI9-BH21'+' '+'BI9-BH19'+' '+'BI9-BH17'+' '+'BI9-BH15'\
          +' '+'BI9-BH13'+' '+'BI9-BH11'+' '+'BI9-BH9'+' '+'BI9-B8'\
          +' '+'BI9-B7.1'+' '+'BI9-B7'+' '+'BI9-BH5'+' '+'BI9-B4'\
          +' '+'BI9-B3'+' '+'BI9-B1'+' '+'G9-BX'+' '+'G10-BX'\
          +' '+'BO10-B1'+' '+'BO10-B3.1'+' '+'BO10-B3.2'+' '+'BO10-B4'\
          +' '+'BO10-BH6'+' '+'BO10-B7'+' '+'BO10-B8'+' '+'BO10-BH10'\
          +' '+'BO10-BH12'+' '+'BO10-BH14'+' '+'BO10-BH16'+' '+'BO10-BH18'\
          +' '+'BO10-BH20'+' '+'BO11-BH20'+' '+'BO11-BH18'+' '+'BO11-BH16'\
          +' '+'BO11-BH14'+' '+'BO11-BH12'+' '+'BO11-BH10'+' '+'BO11-B8'\
          +' '+'BO11-B7'+' '+'BO11-BH6'+' '+'BO11-B4'+' '+'BO11-B3'\
          +' '+'BO11-B1'+' '+'G11-BX'+' '+'G12-BX'+' '+'BI12-B1'\
          +' '+'BI12-B3'+' '+'BI12-B4'+' '+'BI12-BH5'+' '+'BI12-B7'\
          +' '+'BI12-B8'+' '+'BI12-BH9'+' '+'BI12-BH11'+' '+'BI12-BH13'\
          +' '+'BI12-BH15'+' '+'BI12-BH17'+' '+'BI12-BH19'+' '+'BI1-BH21'\
          +' '+'BI1-BH19'+' '+'BI1-BH17'+' '+'BI1-BH15'+' '+'BI1-BH13'\
          +' '+'BI1-BH11'+' '+'BI1-BH9'+' '+'BI1-B8'+' '+'BI1-B7'\
          +' '+'BI1-BH5'+' '+'BI1-B4'+' '+'BI1-B3.2'+' '+'BI1-B3.1'\
          +' '+'BI1-B1'+' '+'G1-BX'+' '+'G2-BX'+' '+'BO2-B1'+' '+'BO2-B3'+\
          ' '+'BO2-B4'+' '+'BO2-BH6'+' '+'BO2-B7'+' '+'BO2-B8'\
          +' '+'BO2-BH10'+' '+'BO2-BH12'+' '+'BO2-BH14'+' '+'BO2-BH16'\
          +' '+'BO2-BH18'+' '+'BO2-BH20'+' '+'BO3-BH20'+' '+'BO3-BH18'\
          +' '+'BO3-BH16'+' '+'BO3-BH14'+' '+'BO3-BH12'+' '+'BO3-BH10'\
          +' '+'BO3-B8'+' '+'BO3-B7.1'+' '+'BO3-B7'+' '+'BO3-BH6'\
          +' '+'BO3-B4'+' '+'BO3-B3'+' '+'BO3-B1'+' '+'G3-BX'+' '+'G4-BX'\
          +' '+'BI4-B1'+' '+'BI4-B3'+' '+'BI4-B4'+' '+'BI4-BH5'\
          +' '+'BI4-B7'+' '+'BI4-B8'+' '+'BI4-BH9'+' '+'BI4-BH11'\
          +' '+'BI4-BH13'+' '+'BI4-BH15'+' '+'BI4-BH17'+' '+'BI4-BH19'\
          +' '+'BI5-BH21'+' '+'BI5-BH19'+' '+'BI5-BH17'+' '+'BI5-BH15'\
          +' '+'BI5-BH13'+' '+'BI5-BH11'+' '+'BI5-BH9'+' '+'BI5-B8'\
          +' '+'BI5-B7'+' '+'BI5-BH5'+' '+'BI5-B4'+' '+'BI5-B3.1'\
          +' '+'BI5-B3'+' '+'BI5-B1'+' '+'G5-BX'
bpmNamesY='G6-BX'+' '+'BO6-B1'+' '+'BO6-B3'+' '+'BO6-B3.1'+' '+'BO6-B4'\
           +' '+'BO6-BV5'+' '+'BO6-B7'+' '+'BO6-B8'+' '+'BO6-BV9'\
           +' '+'BO6-BV11'+' '+'BO6-BV13'+' '+'BO6-BV15'+' '+'BO6-BV17'\
           +' '+'BO6-BV19'+' '+'BO7-BV21'+' '+'BO7-BV19'+' '+'BO7-BV17'\
           +' '+'BO7-BV15'+' '+'BO7-BV13'+' '+'BO7-BV11'+' '+'BO7-BV9'\
           +' '+'BO7-B8'+' '+'BO7-B7'+' '+'BO7-BV5'+' '+'BO7-B4'\
           +' '+'BO7-B3.1'+' '+'BO7-B3'+' '+'BO7-B1'+' '+'G7-BX'\
           +' '+'G8-BX'+' '+'BI8-B1'+' '+'BI8-B3'+' '+'BI8-B3.1'\
           +' '+'BI8-B4'+' '+'BI8-BV6'+' '+'BI8-B7'+' '+'BI8-B8'\
           +' '+'BI8-BV10'+' '+'BI8-BV12'+' '+'BI8-BV14'+' '+'BI8-BV16'\
           +' '+'BI8-BV18'+' '+'BI8-BV20'+' '+'BI9-BV20'+' '+'BI9-BV18'\
           +' '+'BI9-BV16'+' '+'BI9-BV14'+' '+'BI9-BV12'+' '+'BI9-BV10'\
           +' '+'BI9-B8'+' '+'BI9-B7.1'+' '+'BI9-B7'+' '+'BI9-BV6'\
           +' '+'BI9-B4'+' '+'BI9-B3'+' '+'BI9-B1'+' '+'G9-BX'\
           +' '+'G10-BX'+' '+'BO10-B1'+' '+'BO10-B3.1'+' '+'BO10-B3.2'\
           +' '+'BO10-B4'+' '+'BO10-BV5'+' '+'BO10-B7'+' '+'BO10-B8'\
           +' '+'BO10-BV9'+' '+'BO10-BV11'+' '+'BO10-BV13'+' '+'BO10-BV15'\
           +' '+'BO10-BV17'+' '+'BO10-BV19'+' '+'BO11-BV21'\
           +' '+'BO11-BV19'+' '+'BO11-BV17'+' '+'BO11-BV15'\
           +' '+'BO11-BV13'+' '+'BO11-BV11'+' '+'BO11-BV9'+' '+'BO11-B8'\
           +' '+'BO11-B7'+' '+'BO11-BV5'+' '+'BO11-B4'+' '+'BO11-B3'\
           +' '+'BO11-B1'+' '+'G11-BX'+' '+'G12-BX'+' '+'BI12-B1'\
           +' '+'BI12-B3'+' '+'BI12-B4'+' '+'BI12-BV6'+' '+'BI12-B7'\
           +' '+'BI12-B8'+' '+'BI12-BV10'+' '+'BI12-BV12'+' '+'BI12-BV14'\
           +' '+'BI12-BV16'+' '+'BI12-BV18'+' '+'BI12-BV20'+' '+'BI1-BV20'\
           +' '+'BI1-BV18'+' '+'BI1-BV16'+' '+'BI1-BV14'+' '+'BI1-BV12'\
           +' '+'BI1-BV10'+' '+'BI1-B8'+' '+'BI1-B7'+' '+'BI1-BV6'\
           +' '+'BI1-B4'+' '+'BI1-B3.2'+' '+'BI1-B3.1'+' '+'BI1-B1'\
           +' '+'G1-BX'+' '+'G2-BX'+' '+'BO2-B1'+' '+'BO2-B3'+' '+'BO2-B4'\
           +' '+'BO2-BV5'+' '+'BO2-B7'+' '+'BO2-B8'+' '+'BO2-BV9'\
           +' '+'BO2-BV11'+' '+'BO2-BV13'+' '+'BO2-BV15'+' '+'BO2-BV17'\
           +' '+'BO2-BV19'+' '+'BO3-BV21'+' '+'BO3-BV19'+' '+'BO3-BV17'\
           +' '+'BO3-BV15'+' '+'BO3-BV13'+' '+'BO3-BV11'+' '+'BO3-BV9'\
           +' '+'BO3-B8'+' '+'BO3-B7.1'+' '+'BO3-B7'+' '+'BO3-BV5'\
           +' '+'BO3-B4'+' '+'BO3-B3'+' '+'BO3-B1'+' '+'G3-BX'+' '+'G4-BX'\
           +' '+'BI4-B1'+' '+'BI4-B3'+' '+'BI4-B4'+' '+'BI4-BV6'\
           +' '+'BI4-B7'+' '+'BI4-B8'+' '+'BI4-BV10'+' '+'BI4-BV12'\
           +' '+'BI4-BV14'+' '+'BI4-BV16'+' '+'BI4-BV18'+' '+'BI4-BV20'\
           +' '+'BI5-BV20'+' '+'BI5-BV18'+' '+'BI5-BV16'+' '+'BI5-BV14'\
           +' '+'BI5-BV12'+' '+'BI5-BV10'+' '+'BI5-B8'+' '+'BI5-B7'\
           +' '+'BI5-BV6'+' '+'BI5-B4'+' '+'BI5-B3.1'+' '+'BI5-B3'\
           +' '+'BI5-B1'+' '+'G5-BX'

def qStren():
    global initVec
    blah1="cdevCommand "+stone+" get value -names \
    '["+quadNames+"]' | grep ^names > out.txt"
    blah2="cdevCommand "+stone+" get value -names \
    '["+quadNames+"]' | grep ^value >> out.txt"
    blah3="cdevCommand "+stone+" get value -names \
    '["+quadNames+"]' | grep ^trim >> out.txt"
    system(blah1);system(blah2);system(blah3)
    system("sed -e 's/\]//' -e 's/\[//' out.txt > out1.txt")
    system("mv out1.txt out.txt")
    
    initVec=[];initTrim=[]
    x=open('out.txt','r').readlines()
    for i in range(1,len(split(x[1]))):
        #print split(x[0])[i], split(x[1])[i], split(x[2])[i]
        initVec.append(float(split(x[1])[i]))
        initTrim.append(float(split(x[2])[i]))
        
    if (exists('out.txt')):system("rm out.txt")
    return initTrim

def setTrimStren(numQuads,dev):
    list=random.sample(qTrimNames,numQuads);crepe="";qT="";qT0=""
    for j in range(len(list)-1):
        crepe=crepe+list[j]+" ";qT0=qT0+'0.0'+" "
        qT=qT+str( pow(-1,random.randint(1,2))*dev)+" "
    crepe=crepe+list[j+1];qT0=qT0+'0.0'
    qT=qT+str( pow(-1,random.randint(1,2))*dev)
    blah1="cdevCommand "+stone+" set value -names\
    '["+crepe+"]' -trim '["+qT+"]'"
    blah0="cdevCommand "+stone+" set value -names\
    '["+crepe+"]' -trim '["+qT0+"]'"
    
    print blah1;print blah0


def qStrenSet(sign=1.0):
    initTrim=qStren();trimVal=[]
    f=open('/home/cfsb/rcalaga/beta.beat/acd.data/results/qTrim.dat','w')
    x=open('/home/cfsb/rcalaga/beta.beat/acd.data/changeparameters','r').readlines()
    legs=pickle.load(open('/home/cfsb/rcalaga/beta.beat/acd.data/qLengths.pic'))
    
    delta=0
    if len(x)==len(initTrim):
        for j in range(len(x)):
            trim=float(split(x[j])[3])*legs[j]+initTrim[j]
            if abs(trim) < 1e-12:trim=0.0
            
            #trim=float(split(x[j])[3])*legs[j]+initTrim[j]
            trimVal.append(trim);delta=delta+float(split(x[j])[3])*legs[j]
            f.write(split(x[j])[0]+' '+str(initTrim[j]+initVec[j])\
                    +' '+str(trim+initVec[j])+' '+str(initTrim[j])+' '+str(trim)+'\n')
    else: raise ValueError, "Quad Number Mismatch"
    f.close()

    QTRIMS=''
    for k in range(len(trimVal)):     
        QTRIMS=QTRIMS+' '+str(trimVal[k]*sign)
    
    blah1="cdevCommand "+stone+" set value -names\
    '["+quadNames+"]' -trim '["+QTRIMS+"]'"
    print blah1
    print "total KL Change", delta
    #liveStone "set value" -names [bi1-th5 bi1-th3] -trim [.0001 .0002]

def bpmOptics(BNAM,indx):
    #--- get tunes
    crap1='cdevBatch -list "model set default OptiCalc'+r"\n"+\
    stone+' get muX -beamline Blue'+r"\n"+\
    stone+' get muY -beamline Blue" > out.txt';system(crap1)
    x=open('out.txt','r').readlines();tunes=[]
    for i in range(len(x)):
        if len(split(x[i]))>0 and split(x[i])[0] =="value":
            tunes.append(split(x[i])[1])

    #--- get betas & mu's
    blah1='cdevBatch -list "model set default OptiCalc'+r"\n"+\
           stone+' get LatticeFunctions '+\
           '-beamline Blue -names ['+BNAM.lower()+']" > out.txt' 
    system(blah1);system("sed -e 's/\]//' -e 's/\[//' out.txt > out1.txt")
    system("mv out1.txt out.txt")
    
    x=open('out.txt','r').readlines()
    for i in range(1,len(x)):
        if "names" in x[i]: names=x[i]
        if "betax" in x[i]: betx=x[i]
        if "mux" in x[i]:   mux=x[i]
        if "etax" in x[i]:  etax=x[i]
        if "betay" in x[i]: bety=x[i]
        if "muy" in x[i]:   muy=x[i]
        if "scoord" in x[i]:scrd=x[i]
    if indx==0: f=open('twissX.base','w')
    if indx==1: f=open('twissY.base','w')
    f.write('@ Q1               %le    '+tunes[0]+'\n')
    f.write('@ Q2               %le    '+tunes[1]+'\n')
    f.write("* NAME  S  BETX  BETY  MUX   MUY   DX \n")
    f.write("$ %s %le  %le  %le %le  %le  %le \n")
    for j in range(1,len(split(names))):
        f.write('"'+split(names)[j]+'"'+' '+split(scrd)[j]+' '+\
                split(betx)[j]+' '+split(bety)[j]+' '+\
                split(mux)[j]+' '+split(muy)[j]+' '+split(etax)[j]+'\n')
    f.close()
    if (exists('out.txt')):system("rm out.txt")


if __name__ == '__main__':
    print "Calculating Trims..."
    #bpmOptics(bpmNamesX,0);bpmOptics(bpmNamesY,1)
    sign=float(sys.argv[1])
    qStrenSet(sign)

    #--- for trim settings 
    #setTrimStren(3,0.001)


