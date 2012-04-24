#!/afs/cern.ch/eng/sl/lintrack/Python-2.5_32bit/Python-2.5_32bit/bin/python

import  sddsdata as sdds
import re, sys
from metaclass25 import *
from string import *

def Uncrypt(bpm):
  if re.search('/V',bpm):
    plane=1
  elif re.search('/H', bpm):
    plane=0
  else:
    print "No plane found for ", bpm
    plane=-1

  #BPM.20R2.B1/H-bunch110
  bpmname=bpm.split('/')[0]
  bunch=int(bpm.split('bunch')[1])
  return bpmname, plane, bunch


############
#START



filename=sys.argv[1]
print "FILE: ",filename
filenamedir=join(filename.split('/')[0:-1], '/')
if filenamedir=='':
  filenamedir='.'
filenamedir=filenamedir+'/'
onlyfilename=filename.split('/')[-1]


convert=sys.argv[0]  #We assume that twiss.dat is in the same dir as the convert.py
convertdir=join(convert.split('/')[0:-1], '/')
if convertdir=='':
  convertdir='.'
tf=convertdir+'/twiss.dat'
print "TWISS FILE: ", tf 
tw=twiss(tf)


#t=sdds.sddsdata('/afs/cern.ch/user/r/rtomas/lintrack/LHCBPM1-Wed_Sep_10_21-01-36_CEST_2008.sdds.sdds','big')

a=sdds.sddsdata(filename, 'big')

nbOfCapBunches=int(a.data['nbOfCapBunches'])
nbOfCapTurns=int(a.data['nbOfCapTurns'])
print "Turns, Bunches:", nbOfCapTurns,nbOfCapBunches
bpmnames = a.data['bpmNames']
channelnames = a.data['ChannelNames']
#otherBPMnames= [i for i in a.data.keys() if 'BPM' in i]
othernames= [i for i in a.data.keys() if 'BPM' not in i]
#print bpmnames[0],otherBPMnames[0]
#print othernames
#print a.data['BPMSW.1L1.B1/H-bunch111']



bunchfile={}



for bpm in channelnames:
  
    bpmname, plane, bunch = Uncrypt(bpm)
    #print bpm, bpmname, plane, bunch

    try:
      s=tw.S[tw.indx[bpmname]]
    except:
      print bpmname, "Not in model"
      s=-1
             
    if (plane==0 or plane==1) and s > -0.1:
      try:
        bunchfile[bunch]
      except:
        print bunch
        bunchfile[bunch]=open(filenamedir+onlyfilename+str(bunch), "w")
        print >>bunchfile[bunch], "title"


      bunchfile[bunch].write(str(plane)+" "+bpmname+" "+str(s)+" ")
      
      for el in a.data[bpm]:
        bunchfile[bunch].write(str(el)+" ")
      bunchfile[bunch].write("\n")


for i in bunchfile.keys():
  bunchfile[i].close()


#Take the first bunch as default:  
bunch=bunchfile.keys()[0]
os.system('cp '+filenamedir+onlyfilename+' '+filenamedir+onlyfilename+'.sdds')



sys.exit()
print t.data[t.data.keys()[0]]
bpm1=t.data.keys()[1]

print bpm1, t.data[bpm1]
header= t.header
bpmnames = t.data['bpmNames']
print bpmnames
othernames= [i for i in t.data.keys() if 'BPM' not in i]
thevalues = [t.data[i] for i in othernames]
print othernames
print thevalues
#print t.data['bpmNames']
