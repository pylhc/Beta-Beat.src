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
    nturns   = int(d['nbOfCapTurns'].tostring())
    nbunches = int(d['nbOfCapBunches'].tostring())
    print 'I found: ', len(names)
    print nbunches, ' bunches'
    print nturns,   ' turns'
    lbunches = d[plane + 'BunchId'] #has lenght nbunches*nturns*nbpms
    data = d[plane + 'RawData']
    data = data.reshape(nbpms,nturns,nbunches)
    if len(lbunches) / nbpms / nturns != nbunches:
        print 'There is no consistence with header and data hor (nbunches)'
        print len(lbunches_h),nbunches
        print 'First is bunch# ' ,lbunches[0]
        print 'Last  is bunch# ' ,lbunches[len(lbunches)-1]
        bunch = int(lbunches[float(sbunch)]) #will write out table only for selected bunch and plane

    bunch = int(lbunches[float(sbunch)])
    #get bpm position
    vs_twiss, names_twiss = twiss.get_parameter('s'), twiss.get_parameter('name')

    vs= []
    for n in names:
        #print n
        try:
            i = names_twiss.index(n)
            vs.append(float(vs_twiss[i]))
        except IndexError:
            print 'could not find s for ',n,' , setting s to -1000'
            vs.append(-1000.)
    #sorting
    names_sorted=names[argsort(vs)]
    data_sorted =data[argsort(vs)]
    vs_sorted=sorted(vs)

    return [names_sorted,data_sorted,vs_sorted,names,nturns,bunch]

############## => writting to SDDS table  Glenn Vanbavinckhove

####### getting the data
a=sdds.sddsdata(filename, 'big')
a.treat_header()
sddsfile=open(options.PATH+'/'+filename+'.new','w')


###### writing main body
for v in a.par_list:
    txt=v[1].__str__().split()
    txt=txt[0]
    txt=txt.replace("array('","")
    sddsfile.write('# '+v[0]+' '+txt+ '\n')

####### writing main data
#
# plane name pos turn1 turn2 .... turnn

#=> horizontal plane
[names_sorted,data_sorted,vs_sorted,names,nturns,bunch]=getdata('hor')
for name in names_sorted:

   
    bpm = names.tolist().index(name)
    sddsfile.write(str(0)+' '+str(name)+' '+str(twissfile.S[twissfile.indx[name]])+' ')
    for turn in range(nturns):
        sddsfile.write(str(data_sorted[bpm][turn][bunch])+' ')

    sddsfile.write('\n')

    
#=> vertical plane
[names_sorted,data_sorted,vs_sorted,names,nturns,bunch]=getdata('ver')
for name in names_sorted:

   
    bpm = names.tolist().index(name)
    sddsfile.write(str(1)+' '+str(name)+' '+str(twissfile.S[twissfile.indx[name]])+' ')
    for turn in range(nturns):
        sddsfile.write(str(data_sorted[bpm][turn][bunch])+' ')

    sddsfile.write('\n')
    

sddsfile.close()


