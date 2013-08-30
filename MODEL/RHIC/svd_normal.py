#!//afs/cern.ch/eng/sl/lintrack/Python-2.5/bin/python

from os import system
import math, pickle, sys, numpy,re
from string import split, replace
from metaclass import twiss
from operator import mod

def cbpm(s): # some model convention problems
    s=replace(s,"rbpm.","");s=replace(s,"bpm.","")
    s=replace(s,"b-g","g");s=replace(s,"y-g","g")
    s=replace(s,"-bhx","_bx");s=replace(s,"-bvx","_bx")
    if re.match('.*-bh1$',s):s=replace(s,"-bh1","_b1")
    if re.match(".*-bv1$",s):s=replace(s,"-bv1","_b1")
    s=replace(s,"-bh3","_b3");s=replace(s,"-bv3","_b3")
    s=replace(s,"-bh4","_b4");s=replace(s,"-bv4","_b4")
    s=replace(s,"-bh7","_b7");s=replace(s,"-bv7","_b7")
    s=replace(s,"-bh8","_b8");s=replace(s,"-bv8","_b8")
    s=replace(s,"-bh3.1","_b3.1");s=replace(s,"-bv3.1","_b3.1")
    s=replace(s,"-bh3.2","_b3.2");s=replace(s,"-bv3.2","_b3.2")
    s=replace(s,"-bh7.1","_b7.1");s=replace(s,"-bv7.1","_b7.1")
    s=replace(s,"-","_")
    if re.match("^bo10_b3$",s): s=replace(s,"bo10_b3","bo10_b3.1")
    if re.match("^bi1_b3$",s): s=replace(s,"bi1_b3","bi1_b3.1")
    return s.upper()

def writeparams(deltafamilie):
    global variables
    g = open ('changeparameters', 'w');i=0
    for var in variables:
        g.write(var+' ='+ var+'+ ('+str(deltafamilie[i])+');\n')
        i +=1
    g.close()
    return

def model(file):
    mdB=twiss(file);mdlB={}
    for k in range(len(mdB.NAME)):
        mdlB[mdB.NAME[k]]=[mdB.MUX[k],mdB.MUY[k]]
    mdlB['Q1']=mdB.Q1;mdlB['Q2']=mdB.Q2
    return mdlB

def muXY(md,x,y):
    MUX=[];MUY=[];NAMES=[]
    for j in range(len(y.NAME)):
        nm1=cbpm(y.NAME[j]);nm2=cbpm(y.NAME2[j])
        MUY.append(mod(md[nm2][1]-md[nm1][1],1))
    for j in range(len(x.NAME)):
        nm1=cbpm(x.NAME[j]);nm2=cbpm(x.NAME2[j])
        MUX.append(mod(md[nm2][0]-md[nm1][0],1))
    print len(MUX), len(MUY), md['Q1'], md['Q2']
    return MUX,MUY

def twissVec(deltafamilies):
    #print deltafamilies
    writeparams(deltafamilies)
    system('madx < job.iterate.madx > scum');md=model('twiss.dat')
    dir='/afs/cern.ch/user/g/grumolo/scratch1/acd.data/results/'
    x=twiss(dir+'x.out');y=twiss(dir+'y.out')
    #-- read results and make model phadv
    blah=muXY(md,x,y)
    return numpy.concatenate([blah[0],blah[1],[md['Q1']],[md['Q2']]])

def iterate(func, delta1, incr):
    s_matrix=[]
    zerovector=numpy.array(func(delta1))
    pickle.dump(zerovector,open('zerovector','w'),-1)
    #--- calc resp. matrix for each quad circ
    for i in range(0,len(delta1)):
        delta=numpy.array(delta1)
        delta[i]=delta[i]+incr[i]
        print numpy.shape(delta), numpy.shape(incr)
        s_matrix.append(numpy.array((func(delta))-zerovector)/incr[i])
    return numpy.array(s_matrix)


#---- main loop to calc response matrix
def main():
    global variables

    #---open pickle file to read circuits names
    variables=pickle.load(open('qCirc.pic','r'))
    incr=numpy.ones(len(variables),'d')*0.0001
    dfam1=numpy.zeros(len(variables),'d')

    s_matrix=iterate(twissVec,dfam1,incr)
    pickle.dump(s_matrix,open('sensitivity_matrix','w'),-1)

if __name__ == '__main__':
    main()
