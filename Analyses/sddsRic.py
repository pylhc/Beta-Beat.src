#!/afs/cern.ch/eng/sl/lintrack/Python-2.5_32bit/Python-2.5_32bit/bin/python

import  sddsdata as sdds

t=sdds.sddsdata('/afs/cern.ch/user/r/rtomas/lintrack/LHCBPM1-Wed_Sep_10_21-01-36_CEST_2008.sdds.sdds','big')


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
