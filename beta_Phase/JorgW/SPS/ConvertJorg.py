

import gzip
from optparse import OptionParser
import os, sys
from re import match
from YASPmetaclassSPS import *

parser = OptionParser()
parser.add_option("-a", "--accel",
                help="What accelerator: LHCB1 LHCB2 SPS RHIC",
                metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-p", "--path",
                help="Path to data",
                metavar="PATH", default="./",dest="path")
parser.add_option("-o", "--out",
                help="Output path",
                metavar="OUT", default="./",dest="output")
parser.add_option("-x", "--qx",
                help="Horizontal tune",
                metavar="Q1", default="0",dest="q1")
parser.add_option("-y", "--qy",
                help="Vertical tune",
                metavar="Q2", default="0",dest="q2")
parser.add_option("-d", "--dictionary",
                help="Dictionary for BPM names between YASP and MODEL",
                metavar="DICT", default="0",dest="dict")

parser.add_option("-m", "--model",
                help="Twiss model to run Akio's code",
                metavar="Twiss", default="0",dest="twiss")


(options, args) = parser.parse_args()

p=sys.argv[0]
execpath =p[:-len(p.split("/")[-1])]
print "Path to the executable", execpath
Akiocode=execpath+"../../SAD/BetaAnalysis.ro.sad"

path=options.path+'/'
output=options.output+'/'
q1=options.q1
q2=options.q2

if q1=="0":
  print "Warning: Horizontal tune not given in input"
  q1=0.13
if q2=="0":
  print "Warning: Vertical tune not given in input"
  q2=0.18
if options.dict=="0":
  YASPdictionary={}
else:
  execfile(options.dict)


files=os.listdir(path)
files=filter(lambda x: "RM" in x, files) # Assume all YASP files contain "RM"
files=filter(lambda x: "_lin" not in x, files) # Never use already existing *_lin[xy] files!


a=[]
samecorr={}
onlystr=0
for filename in files:
    a.append(YASPtwiss(path+filename, YASPdictionary))
    
    print filename, a[-1].exciter, a[-1].plane, a[-1].str 
    if abs(float(a[-1].str)) > onlystr:  #Only use cases with |str| > onlystr
      try:
        samecorr[a[-1].exciter].append(len(a)-1)
      except:
        samecorr[a[-1].exciter]=[]
        samecorr[a[-1].exciter].append(len(a)-1)
    
print samecorr


for corrs in samecorr:
    jn=samecorr[corrs]
    
    
    foutx=open(output+corrs+'_linx', "w")
    fouty=open(output+corrs+'_liny', "w")
    print "Writing to ", output+corrs
    print >> foutx, "@ DIR %01s \""+a[jn[0]].plane+"\""
    print >> foutx, "@ EXCITER %s \""+corrs+"\""
    print >> foutx, "@ Q1 %le ", q1
    print >> foutx, "* NAME   CO"
    print >> foutx, "$ %s   %le"

    print >> fouty, "@ DIR %01s \""+a[jn[0]].plane+"\""
    print >> fouty, "@ EXCITER %s \""+corrs+"\""
    print >> fouty, "@ Q2 %le ", q2
    print >> fouty, "* NAME   CO"
    print >> fouty, "$ %s   %le"

    for name in a[jn[0]].HNAME:
      orb=0
      problem=0
      for ii in jn:
        try:
          iin=a[ii].Hindx[name]
          orb = orb +  a[ii].HCO[iin] * a[ii].str
        except:
          print "problem with ", name
          problem=1
      if problem==0:
        print >> foutx, "\""+name+"\"", orb
    
    for name in a[jn[0]].VNAME:
      orb=0
      problem=0
      for ii in jn:
        try:
          iin=a[ii].Vindx[name]
          orb = orb +  a[ii].VCO[iin] * a[ii].str
        except:
          print "problem with ", name
      if problem==0:
        print >> fouty, "\""+name+"\"", orb

    
    foutx.close()
    fouty.close()
    

command=Akiocode+" -o "+output+"/beta.dat "+options.twiss+" "+output+"/*_lin[xy]"
print "Running Akio's code with command:", command



stat=os.system(command)
print "STATUS ", stat




from re import match
from metaclass import *

a=twiss(output+'/beta.dat')
b=twiss(options.twiss)

fx=open(output+"/getphasex.out", "w")
fy=open(output+"/getphasey.out", "w")

print >> fx, "@ Q1 %le ",q1
print >> fy, "@ Q2 %le ",q2

print >> fx,"* NAME NAME2 POS1 POS2 PHASE RMS PHASEMDL"
print >> fx,"$ %s   %s    %le  %le  %le  %le  %le"

print >> fy,"* NAME NAME2 POS1 POS2 PHASE RMS PHASEMDL"
print >> fy,"$ %s   %s    %le  %le  %le   %le %le"


indx=[]
indy=[]

for i in range(len(a.NAME)):
    if a.MUX[i] >  -1:
        indx.append(i)
    if a.MUY[i] >  -1:
        indy.append(i)


for ix in range(1,len(indx)):
    ix1=indx[ix-1]
    ix2=indx[ix]
    indm1=b.indx[a.NAME[ix1]]
    indm2=b.indx[a.NAME[ix2]]
    print >> fx, "\""+a.NAME[ix1]+"\"" , "\""+a.NAME[ix2]+"\"", a.S[ix1], a.S[ix2], a.MUX[ix2]-a.MUX[ix1], 0, b.MUX[indm2]-b.MUX[indm1]



for ix in range(1,len(indy)):
    ix1=indy[ix-1]
    ix2=indy[ix]
    indm1=b.indx[a.NAME[ix1]]
    indm2=b.indx[a.NAME[ix2]]
    print >> fy, "\""+a.NAME[ix1]+"\"" , "\""+a.NAME[ix2]+"\"", a.S[ix1], a.S[ix2], a.MUY[ix2]-a.MUY[ix1], 0, b.MUY[indm2]-b.MUY[indm1]

