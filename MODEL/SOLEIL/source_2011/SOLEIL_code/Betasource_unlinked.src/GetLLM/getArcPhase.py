from metaclass import *
from optparse import OptionParser
import os,re,math
from numpy.linalg import *

parser = OptionParser()


parser.add_option("-p", "--path",
        help="Path to measurement",
        metavar="PATH", default="0",dest="path")

parser.add_option("-o", "--out",
        help="output path",
        metavar="OUT", default="./",dest="out")


parser.add_option("-b", "--bpms",
        help="bpm file to compute phase differences",
        metavar="BPMS", default="/afs/cern.ch/user/r/rtomas/lintrack/Beta-Beat.src/GetLLM/bpms.out",dest="bpms")


parser.add_option("-l", "--label",
        help="Label for output files",
        metavar="LABEL", default="",dest="label")




(options, args) = parser.parse_args()


os.system("mkdir -p "+options.out)

label=options.label
outx=open(options.out+"/arcphasediffx_"+label+".dat","w")
outy=open(options.out+"/arcphasediffy_"+label+".dat","w")

out1x=open(options.out+"/arcphasex_"+label+".dat","w")
out1y=open(options.out+"/arcphasey_"+label+".dat","w")

out2x=open(options.out+"/phasediffx_"+label+".dat","w")
out2y=open(options.out+"/phasediffy_"+label+".dat","w")


try:
  tx=twiss(options.path+'/getphasetotx_free.out')
  ty=twiss(options.path+'/getphasetoty_free.out')
except:
  tx=twiss(options.path+'/getphasetotx.out')
  ty=twiss(options.path+'/getphasetoty.out')



if 'B1' in tx.NAME[0]:
  beam='b1'
  b='B1'
else:
  beam='b2'
  b='B2'

print "This is ", beam


def aveArcx(IP, side , t):
  #m1=re.compile('\.1[0-3]'+side+IP)
  m2=re.compile('\.[5-9]'+side+IP)
  m1=re.compile('\.1[0-3]NO_BPM_WILL_MATCH'+side+IP)  # Very ugly way to avoid using m1
  #m2=re.compile('\.[4-8]'+side+IP)
  ave=0
  ave2=0
  count=0
  s=0
  for i in range(len(t.NAME)):
    if m1.search(t.NAME[i]) or m2.search(t.NAME[i]):

      ph=(t.PHASEX[i]-t.PHXMDL[i])
      if ph<-0.5: ph=ph+1
      if ph>0.5: ph=ph-1
      s=s+t.S[i]
      ave=ave+ph
      ave2=ave2+ph**2
      count=count+1
      
  if count>0:
      ave=ave/count
      ave2=ave2/count
      return count, ave, math.sqrt(ave2-ave**2+1e-16), s/count
  else:
      return 0,0,0,0


def aveArcy(IP, side , t):
    #m1=re.compile('\.1[0-3]'+side+IP)
    m2=re.compile('\.[5-9]'+side+IP)
    m1=re.compile('\.1[0-3]NO_BPM_WILL_MATCH'+side+IP)  # Very ugly way to avoid using m1
    #m2=re.compile('\.[6-9]'+side+IP)
    ave=0
    ave2=0
    count=0
    s=0
    for i in range(len(t.NAME)):
        if m1.search(t.NAME[i]) or m2.search(t.NAME[i]):
            ph=(t.PHASEY[i]-t.PHYMDL[i])
            if ph<-0.5: ph=ph+1
            if ph>0.5: ph=ph-1
            s=s+t.S[i]
            ave=ave+ph
            ave2=ave2+ph**2
            count=count+1
    if count>0:
      ave=ave/count
      ave2=ave2/count
      return count, ave, math.sqrt(ave2-ave**2+1e-16),s/count
    else:
      return 0,0,0,0



resx={}
resy={}
for IP in ['1','2','3','4','5','6','7','8']:
  for side in ['L','R']:
    
    resx[side+IP]=aveArcx(IP,side,tx)
    print >>out1x, side+IP,  resx[side+IP][3],  resx[side+IP][1], resx[side+IP][2], resx[side+IP][0]
    resy[side+IP]=aveArcy(IP,side,ty)
    print >>out1y, side+IP, resy[side+IP][3],  resy[side+IP][1], resy[side+IP][2], resy[side+IP][0]


out1x.close()
out1y.close()

bpms=twiss(options.bpms)

def phasediffx(i,j):
  ph=(tx.PHASEX[i]-tx.PHASEX[j])
  if ph<-0.5: ph=ph+1
  if ph>0.5: ph=ph-1
  phmdl=(tx.PHXMDL[i]-tx.PHXMDL[j])
  if phmdl<-0.5: phmdl=phmdl+1
  if phmdl>0.5: phmdl=phmdl-1
  std = sqrt(tx.STDPHX[i]**2 + tx.STDPHX[j]**2)
  return tx.NAME[i],tx.NAME[j],ph,std,phmdl

def phasediffy(i,j):
  ph=(ty.PHASEY[i]-ty.PHASEY[j])
  if ph<-0.5: ph=ph+1
  if ph>0.5: ph=ph-1
  phmdl=(ty.PHYMDL[i]-ty.PHYMDL[j])
  if phmdl<-0.5: phmdl=phmdl+1
  if phmdl>0.5: phmdl=phmdl-1
  std = sqrt(ty.STDPHY[i]**2 + ty.STDPHY[j]**2)
  return ty.NAME[i],ty.NAME[j],ph,std,phmdl

resIPx={}
resIPy={}
for i in range(0,len(bpms.NAME1)):
  resIPx[i]=phasediffx(tx.indx[bpms.NAME1[i]+'.'+b],tx.indx[bpms.NAME2[i]+'.'+b])
  print >>out2x, resIPx[i][0],resIPx[i][1],resIPx[i][2],resIPx[i][3],resIPx[i][4],  resIPx[i][4]-resIPx[i][2], resIPx[i][3]

for i in range(0,len(bpms.NAME1)):
  resIPy[i]=phasediffy(ty.indx[bpms.NAME1[i]+'.'+b],ty.indx[bpms.NAME2[i]+'.'+b])
  print >>out2y, resIPy[i][0],resIPy[i][1],resIPy[i][2],resIPy[i][3],resIPy[i][4],  resIPy[i][4]-resIPy[i][2], resIPy[i][3]

out2x.close()
out2y.close()



print >>outx, "* ARC DPHASE ERR"
print >>outx, "$  %s  %le %le"

print >>outy, "* ARC DPHASE ERR"
print >>outy, "$  %s  %le %le"

#resp matrix
# phix = 0.003262841/0.0005 *kqtd.a67b2 + 0.017476/0.0005 *kqtf.a67b2 
# phiy =-0.01805/0.0005     *kqtd.a67b2 + -0.00339/0.0005 *kqtf.a67b2


if beam=='b1':

  R=[[0.017476/0.0005,0.003262841/0.0005],[-0.00339/0.0005,-0.01805/0.0005]]

if beam=='b2':

  R=[[0.0177/0.0005,0.0033262/0.0005],[-0.00333/0.0005,-0.01785/0.0005]]




Ri=inv(array(R))
print Ri
# Very silly thing to avoid problems in numpy version of dev server!!!
Ri=[Ri[0][0],Ri[0][1] ], [Ri[1][0],Ri[1][1] ]
print Ri

totx=0
toty=0
qf={}
qd={}
for arc in [['1','2'],['2','3'],['3','4'],['4','5'],['5','6'],['6','7'],['7','8'],['8','1']]:
  a=arc[0]+arc[1]
  L='L'+arc[1]
  R='R'+arc[0]
  x=resx[L][1]-resx[R][1]
  y=resy[L][1]-resy[R][1]
  totx=totx+x
  toty=toty+y
  [qf[a],qd[a]]=dot(Ri,[x,y])
  print >>outx, a, x, math.sqrt(resx[L][2]**2+resx[R][2]**2)
  print >>outy, a, y, math.sqrt(resy[L][2]**2+resy[R][2]**2)


print >> outx, "@ TOT %le ", totx
print >> outy, "@ TOT %le ", toty






[qf['all'],qd['all']]=dot(Ri,array([totx,toty]))

print qf,qd

for arc in [['2','3'],['3','4'],['6','7'],['7','8']]:
    a=arc[0]+arc[1]
    qf[a]=qf[a]-qf['all']/4
    qd[a]=qd[a]-qd['all']/4
    



f=open("changeparameters.tfs","w")
fm=open("changeparameters.madx","w")

print >>f, "* NAME DELTA"
print >>f, "$  %s  %le"
for arc in [['1','2'],['2','3'],['3','4'],['4','5'],['5','6'],['6','7'],['7',   '8'],['8','1']]:
    a=arc[0]+arc[1]
    print >>f,'kqtf.a'+a+beam, -qf[a]
    print >>f,'kqtd.a'+a+beam, -qd[a]
    print >>fm,'kqtf.a'+a+beam,'=', 'kqtf.a'+a+beam,'-(',qf[a],');'
    print >>fm,'kqtd.a'+a+beam,'=', 'kqtd.a'+a+beam,'-(',qd[a],');'

    





f.close()
outx.close()
outy.close()




