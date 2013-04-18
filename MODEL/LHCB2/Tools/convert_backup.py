#F.R. 25-Aug-09

# => changed bu Glenn Vanbavinckhove to make it compatible

import  sddsdataFED as sdds
import io,sys
from numpy import *
from metaclass import twiss as twiss2
from optparse import OptionParser

# example : /afs/cern.ch/eng/sl/lintrack/Python-2.5_32bit/Python-2.5_32bit/bin/python lhcconverter.py -f LHCBPM2-Mon_Jun_29_16-26-00_CEST_2009.sdds -b 3 -m /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB2/twiss.dat -p ./

parser = OptionParser()
parser.add_option("-f", "--file", 
		 help="Filename that has to be converted",
		 metavar="FILEN", default="null",dest="FILEN")
parser.add_option("-b", "--bunch", 
		 help="Which bunch to use",
		 metavar="BUNCH", default="0",dest="BUNCH")
parser.add_option("-m", "--model", 
		 help="Which bunch to use",
		 metavar="MODEL", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB1/nominal.opt",dest="MODEL")
parser.add_option("-p", "--poutput", 
		 help="Output path",
		 metavar="PATH", default="./",dest="PATH")


            


(options, args) = parser.parse_args()

##########

filename=options.FILEN
twiss = io.TFSReader(options.MODEL)
twissfile=twiss2(options.MODEL)
sbunch = options.BUNCH #this is index of list of bunches in lbunches (-->bunch id array) 

#ex: lbunches=['3','4','5','6'], sbunch=3 --> bunch_id=6


##############################
#load sdds data
def getdata(plane):
    
   
    print 'data list has ',len(a.data), ' elements'
    d = a.data[0]
    names    = d['bpmNames']
    nbpms    = len(names)
    print "Number of bpms ",nbpms
    #print "the index of the bpms is ",names.tolist().index('BPM.29R4.B1')
    #print d['nbOfCapTurns']
    #nturns   = int(d['nbOfCapTurns'].tostring())
    #nbunches = int(d['nbOfCapBunches'].tostring())
    nturns=int(nturnss)
    nbunches=int(nbunchess)
    print 'I found: ', len(names)
    print nbunches, ' bunches'
    print nturns,   ' turns'
    lbunches = d[plane + 'BunchId'] #has lenght nbunches*nturns*nbpms
    data = d[plane + 'PositionsConcentratedAndSorted']
    #print data
    #data = data.reshape(nbpms,nturns,nbunches)
    data = data.reshape(nbpms,nbunches,nturns)
    #print data
    print 'First is bunch# ' ,lbunches[0]
    print 'Last  is bunch# ' ,lbunches[len(lbunches)-1]
    #for i in lbunches:
     #   print i
    if len(lbunches) / nbpms / nturns != nbunches:
        print 'There is no consistence with header and data for (nbunches)'
        print len(lbunches_h),nbunches
        print 'First is bunch# ' ,lbunches[0]
        print 'Last  is bunch# ' ,lbunches[len(lbunches)-1]
        #bunch = int(lbunches[float(sbunch)]) #will write out table only for selected bunch and plane

    #print lbunches[len(lbunches)-1]
    bunch='none'
    for bunchh in range(0,int(nbunchess)):

        #print lbunches[bunchh],sbunch
        if int(lbunches[bunchh])==int(sbunch):

            bunch=bunchh
            print 'bunch found at position ',bunch,sbunch,lbunches[bunchh]
            break            

    if bunch=='none':
        print "no bunch found"
        sys.exit()
  
        
    #print "Looking in bunch ", bunch 
    #get bpm position
    vs_twiss, names_twiss = twiss.get_parameter('s'), twiss.get_parameter('name')

    vs= []
    for n in names:
        #print n
        try:
            i = names_twiss.index(n)
            vs.append(float(vs_twiss[i]))
            #print vs_twiss[i]
        except IndexError:
            print 'could not find s for ',n,' , setting s to -1000'
            vs.append(-1000.)
    #sorting
    names_sorted=names[argsort(vs)]
    data_sorted =data[argsort(vs)]
    #for i in range(0,10):
        #print data[94]
    vs_sorted=sorted(vs)

    return [names_sorted,data_sorted,vs_sorted,names,nturns,bunch]


############## => writting to SDDS table  Glenn Vanbavinckhove

####### getting the data
a=sdds.sddsdata(filename, 'big')
a.treat_header()
sddsfile=open(filename+'.sdds.new','w')
nturns=''
nbunches=''

#print "Loading data in bunch ",bunch

###### writing main body
for v in a.par_list:
    #print v
    txt=v[1].__str__().split()
    txt=txt[0]
    txt=txt.replace("array('","")
    if v[0]=='nbOfCapTurns':
	nturnss=txt
    elif v[0]=='nbOfCapBunches':
	nbunchess=txt
        #print "bunches "+txt
    #print txt
    sddsfile.write('# '+v[0]+' '+txt+ '\n')

#sddsfile=open(filename+'.sdds.new','w')
####### writing main data
#
# plane name pos turn1 turn2 .... turnn

#=> horizontal plane
[names_sorted,data_sorted,vs_sorted,names,nturns,bunch]=getdata('hor')
for name in names_sorted:

   
    bpm = names.tolist().index(name)
    sddsfile.write(str(0)+' '+str(name)+' '+str(twissfile.S[twissfile.indx[name]])+' ')
    for turn in range(nturns):
        sddsfile.write(str(data_sorted[bpm][int(bunch)][turn])+' ')

    sddsfile.write('\n')

    
#=> vertical plane
[names_sorted,data_sorted,vs_sorted,names,nturns,bunch]=getdata('ver')

for name in names_sorted:

   
    bpm = names.tolist().index(name)
    sddsfile.write(str(1)+' '+str(name)+' '+str(twissfile.S[twissfile.indx[name]])+' ')
    for turn in range(nturns):
        sddsfile.write(str(data_sorted[bpm][int(bunch)][turn])+' ')

    sddsfile.write('\n')
    

sddsfile.close()


