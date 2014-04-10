import sys

import json

#execfile('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/LHCB1/AllLists.py') 
#vco,vb1=getListsByIR()
vco,vb1=json.load(file('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/LHCB1/AllLists.json','r'))['getListsByIR']


#execfile('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/LHCB2/AllLists.py') 
#vco,vb2=getListsByIR()
vco,vb2=json.load(file('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/LHCB2/AllLists.json','r'))['getListsByIR']



IPN=1
if len(sys.argv) > 1:
 IPN=int(sys.argv[1])

if IPN>8 or IPN<1:
 IPN=1

#IPNO='\''+str(IPN)+'\''
IPNO=str(IPN)
#print IPN, IPNO


fouta= open("applycorrection.seqx",'w')

foutc= open("variablesc.seqx",'w')
foutb1= open("variablesb1.seqx",'w')
foutb2= open("variablesb2.seqx",'w')
fout2 = open("svariables.seqx",'w')
fout3 = open("dvariables.seqx",'w')

fout = open("genchangpars.seqx",'w')
fout.write('select,flag=save, clear;')


v = vb1[IPNO];
print '\nBeam 1\n'
fout.write('!B1\n'); 
for q in v:
  print q
  foutb1.write('   vary, name=d'+q+', step:=1e-4;\n');
  fout2.write(' '+q+'_0 = '+q+';\n');
  fout3.write(' '+q+' := '+q+'_0 + d'+q+';\n');
  fout.write('select,flag=save,pattern=\"d'+q+'\";\n')
  
  fouta.write(q+' = '+q+'_0 + d' +  q+';\n');
################################################################


v = vb2[IPNO];
print '\nBeam 2\n'
fout.write('\n!B2\n'); 
for q in v:
  print q
  foutb2.write('   vary, name=d'+q+', step:=1e-4;\n');
  fout2.write(' '+q+'_0 = '+q+';\n');
  fout3.write(' '+q+' := '+q+'_0 + d'+q+';\n');
  fout.write('select,flag=save,pattern=\"d'+q+'\";\n')
  fouta.write(q+' = '+q+'_0 + d' +  q+';\n');

################################################################

v = vco[IPNO];
print '\nBeam 1 and Beam 2\n'
fout.write('\n!B1 and B2\n'); 
for q in v:
  print q
  foutc.write('   vary, name=d'+q+', step:=1e-4;\n');
  fout2.write(' '+q+'_0 = '+q+';\n');
  fout3.write(' '+q+' := '+q+'_0 + d'+q+';\n');
  fout.write('select,flag=save,pattern=\"d'+q+'\";\n')
  fouta.write(q+' = '+q+'_0 + d' +  q+';\n');
  
foutc.close();
foutb1.close();
foutb2.close();
fout2.close();  
fout3.close();

fout.write('\n save, file=\"changeparameters.madx\";\n');
fout.close()
