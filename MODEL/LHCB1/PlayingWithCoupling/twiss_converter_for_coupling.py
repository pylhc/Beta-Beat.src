from metaclass import *


twissin=sys.argv[1]
a=twiss(twissin)
a.Cmatrix()

f=open('twiss.couple.rm.dat','w')
g=open(twissin,'r')

for line in g:
    sl=line.split(" ")
    if line[0]=="@":
        f.write(line)

f.write('* NAME   S      BETX   BETY   MUX    MUY    X      Y      DX     ALFX   ALFY   F1001  F1001R F1001I F1010  F1010R F1010I\n')
f.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le\n')

for i in range(len(a.NAME)):
    f.write('"'+a.NAME[i]+'" '+str(a.S[i])+' '+str(a.BETX[i])+' '+str(a.BETY[i])+' '+str(a.MUX[i])+' '+str(a.MUY[i])+' '+str(a.X[i])+' '+str(a.Y[i])+' '+str(a.DX[i])+' '+str(a.ALFX[i])+' '+str(a.ALFY[i])+' '+str(abs(a.f1001[i]))+' '+str(a.f1001[i].real)+' '+str(a.f1001[i].imag)+' '+str(abs(a.f1010[i]))+' '+str(a.f1010[i].real)+' '+str(a.f1010[i].imag)+'\n')


f.close()
