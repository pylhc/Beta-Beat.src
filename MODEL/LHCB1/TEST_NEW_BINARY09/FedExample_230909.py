#F.R. 25-Aug-09
#to do: 
#-chack with real orbit!
#-maybe adjust space with .rjust,.ljust for better layout of tfs output
#-once there is data, verify we have enough digits on position data printed out in tfs
#-could setup:
#---- read input info from argv
#---- read a list of bunches insteda of a single one for run
#---- could change input sbunch to 
#     'bunch_id' (now it's just index of the list of available bunches)

import  sddsdataFED as sdds
import io
from numpy import *

#input, should be set up with argv from the command line?

filename="LHCBPM2-Mon_Jun_29_16-26-00_CEST_2009.sdds"
twiss = io.TFSReader('../../LHCB2/twiss.dat')

#select the plane
plane='hor'
#plane='ver'

#select bunches for which you want to print out the orbit
sbunch = 3 #this is index of list of bunches in lbunches (-->bunch id array) 

#ex: lbunches=['3','4','5','6'], sbunch=3 --> bunch_id=6


##############################
#load sdds data

a=sdds.sddsdata(filename, 'big')


print 'data list has ',len(a.data), ' elements'

#assuming that there is 1 single dictionary:
d = a.data[0]

names    = d['bpmNames']
nbpms    = len(names)
nturns   = int(d['nbOfCapTurns'].tostring())
nbunches = int(d['nbOfCapBunches'].tostring())


print 'I found:\n ', len(names), ' BPMs: ',names
print nbunches, '\t bunches'
print nturns,   '\t turns'

### data


lbunches = d[plane + 'BunchId'] #has lenght nbunches*nturns*nbpms

data = d[plane + 'RawData']
#data = d[plane + 'PositionsConcentratedAndSorted'] #need to use this???
data = data.reshape(nbpms,nturns,nbunches)

#lbunches.sort()

if len(lbunches) / nbpms / nturns != nbunches:
    print 'There is no consistence with header and data hor (nbunches)'
    print len(lbunches_h),nbunches

print 'First is bunch# ' ,lbunches[0]
print 'Last  is bunch# ' ,lbunches[len(lbunches)-1]



bunch = int(lbunches[sbunch]) #will write out table only for selected bunch and plane

#########initialize tfs output
t=io.TFSWriter()

#treat sdds header parameters for creating TFS header
a.treat_header()

for v in a.par_list:
    vv=[]
    for val in v:
        vv.append(val.__str__().replace('\n',''))
    t.add_scalar_parameter_value(vv[0].ljust(18),vv[1].ljust(8),vv[2].rjust(1))

###########
outname = filename.replace('.sdds','_BUNCH_'+str(bunch)+'_'+plane+'.tfs')

#get bpm position
vs_twiss, names_twiss = twiss.get_parameter('s'), twiss.get_parameter('name')

vs= []
for n in names:
    try:
        i = names_twiss.index(n)
        vs.append(float(vs_twiss[i]))
    except IndexError:
        print 'could not find s for ',n,' , setting s to -1000'
        vs.append(-1000.)
   

#create variable names for TFS table columns
col_names, col_types = [],[]

col_names.append('BPM')
col_types.append('%s')

col_names.append('s')
col_types.append('%le')

cwidth = len('H-orbit-Turn#1000')
for n in range(nturns):
    col_names.append('H-orbit-Turn#'+ str(n+1).rjust(cwidth))
    col_types.append('%le'.rjust(cwidth))

t.set_parameters(col_names, col_types)


#sorting
names_sorted=names[argsort(vs)]
data_sorted =data[argsort(vs)]

vs_sorted=sorted(vs)
#write out twiss table (bpmname,s,turn 1...turn n)
for s,name in zip(vs_sorted,names_sorted):
    record = []
    record.append(name)
    record.append(str(round(float(s),3)).rjust(10))
    bpm = names.tolist().index(name)
    for turn in range(nturns):
        record.append(str(data_sorted[bpm][turn][bunch]))
    t.add_record(record)
####
t.write(outname)



