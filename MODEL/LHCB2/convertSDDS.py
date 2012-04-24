
#
# example:
#
# python convertSDDS.py /nfs/cs-ccr-nfs4/slopsdata/2008/OP_DATA/LHCBPM-SDDS/LHCBPM2-Thu_Sep_04_22-09-52_CEST_2008.sdds.gz


from sddsRam import *
import re, sys
from metaclass import *
from string import *
filename=sys.argv[1]
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
tw=twiss(tf)


a=sddsRead(filename)


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


nbOfCapBunches=int(a.data['nbOfCapBunches'])
nbOfCapTurns=int(a.data['nbOfCapTurns'])

print "Turns, Bunches:", nbOfCapTurns,nbOfCapBunches

bunchfile={}


for bpm in a.data.keys():
  if re.search('BPM', bpm):
    bpmname, plane, bunch = Uncrypt(bpm)
    #print bpm, bpmname, plane, bunch

    try:
      s=tw.S[tw.indx[bpmname]]
    except:
      print bpmname, "Not in model"
             
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
