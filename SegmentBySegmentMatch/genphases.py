import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD") 
from metaclass import *
import operator


IPN=1
if len(sys.argv) > 1:
 IPN=int(sys.argv[1])

if IPN>8 or IPN<1:
 IPN=1

IPNO=str(IPN)
print IPN

#datab1 = twiss('Beam1/sbs/twiss_IP'+IPNO+'.dat')

#              'Beam1/sbs/sbsphasext_IP1.out')
datab1 = twiss('Beam1/sbs/sbsphasext_IP'+IPNO+'.out')
#print len(datab1.NAME)
#print datab1.NAME[0]

datab2 = twiss('Beam2/sbs/sbsphasext_IP'+IPNO+'.out')
#print len(datab2.NAME)
#print datab2.NAME[0]

import numpy
s1 = datab1.MODEL_S
indsb1 = numpy.argsort(s1)
print indsb1

s2 = datab2.MODEL_S
indsb2 = numpy.argsort(s2)
print indsb2

print 'Start H B1: ',datab1.NAME[indsb1[0]]
print 'Start H B2: ',datab2.NAME[indsb2[0]]

#print s1[inds[0]]
#print s1[inds[1]]
#print s1[inds[2]]
#sys.exit();


fout1 = open("phases.seqx",'w')
fout0b1 = open("phases0b1.seqx",'w')
fout0b2 = open("phases0b2.seqx",'w')



print datab1.NAME[0], datab1.NAME[len(datab1.NAME)-1]
print datab2.NAME[0], datab2.NAME[len(datab2.NAME)-1]

##########################################################

fout1.write('\n!!!! BEAM 1 H !!!!!\n\n'); 

fout0b1.write('\n!!!! BEAM 1 H !!!!!\n\n'); 

for name in datab1.NAME:
# print name

 fout0b1.write('mux0'+name+' = ');
 fout0b1.write('table(twiss, '+name+', mux) - ');
 fout0b1.write('table(twiss, '+datab1.NAME[indsb1[0]]+', mux);\n');

 fout1.write('mux'+name+' := ');
 fout1.write('table(twiss, '+name+', mux) - ');
 fout1.write('table(twiss, '+datab1.NAME[indsb1[0]]+', mux);\n');


fout1.write('\n');
for name in datab1.NAME:
 fout1.write('dmux'+name+' := ');
 fout1.write('mux'+name+' - ' 'mux0'+name+';\n'  );


##########################################################

fout1.write('\n!!!! BEAM 2 H !!!!!\n\n'); 

fout0b2.write('\n!!!! BEAM 2 H !!!!!\n\n'); 

for idx in indsb2:
# print idx
 name = datab2.NAME[idx]
# print name

 fout0b2.write('mux0'+name+' = ');
 fout0b2.write('table(twiss, '+name+', mux) - ');
 fout0b2.write('table(twiss, '+datab2.NAME[indsb2[0]]+', mux);\n');

 fout1.write('mux'+name+' := ');
 fout1.write('table(twiss, '+name+', mux) - ');
 fout1.write('table(twiss, '+datab2.NAME[indsb2[0]]+', mux);\n');
 

###
###


fout1.write('\n');
for name in datab2.NAME:
 fout1.write('dmux'+name+' := ');
 fout1.write('mux'+name+' - ' 'mux0'+name+';\n'  );
  

##########################################################

datab1 = twiss('Beam1/sbs/sbsphaseyt_IP'+IPNO+'.out')
datab2 = twiss('Beam2/sbs/sbsphaseyt_IP'+IPNO+'.out') #Beam2/sbs/sbsphaseyt_IP1.out

s1 = datab1.MODEL_S
indsb1 = numpy.argsort(s1)
print indsb1

s2 = datab2.MODEL_S
indsb2 = numpy.argsort(s2)
print indsb2

print 'Start V B1: ',datab1.NAME[indsb1[0]]
print 'Start V B2: ',datab2.NAME[indsb2[0]]

fout1.write('\n!!!! BEAM 1 V !!!!!\n\n'); 

fout0b1.write('\n!!!! BEAM 1 V !!!!!\n\n'); 

for name in datab1.NAME:
# print name

 fout0b1.write('muy0'+name+' = ');
 fout0b1.write('table(twiss, '+name+', muy) - ');
 fout0b1.write('table(twiss, '+datab1.NAME[indsb1[0]]+', muy);\n');

 fout1.write('muy'+name+' := ');
 fout1.write('table(twiss, '+name+', muy) - ');
 fout1.write('table(twiss, '+datab1.NAME[indsb1[0]]+', muy);\n');


fout1.write('\n');
for name in datab1.NAME:
 fout1.write('dmuy'+name+' := ');
 fout1.write('muy'+name+' - ' 'muy0'+name+';\n'  );


##########################################################

fout1.write('\n!!!! BEAM 2 V !!!!!\n\n'); 

fout0b2.write('\n!!!! BEAM 2 V !!!!!\n\n'); 

for name in datab2.NAME:
# print idx
 print name

 fout0b2.write('muy0'+name+' = ');
 fout0b2.write('table(twiss, '+name+', muy) - ');
 fout0b2.write('table(twiss, '+datab2.NAME[indsb2[0]]+', muy);\n');

 fout1.write('muy'+name+' := ');
 fout1.write('table(twiss, '+name+', muy) - ');
 fout1.write('table(twiss, '+datab2.NAME[indsb2[0]]+', muy);\n');



fout1.write('\n');
for name in datab2.NAME:
 fout1.write('dmuy'+name+' := ');
 fout1.write('muy'+name+' - ' 'muy0'+name+';\n'  );
  

##########################################################


fout0b1.close();
fout0b2.close();
fout1.close();



#sigmaStart=[0.000227,9.306e-5, 5.02e-4, 4.21e-5,0.00666,0.002]
#map = Map(4,'theTransferMap.map');
