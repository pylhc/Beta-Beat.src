from pylab import *
from Timber import *
from optparse import OptionParser
import numpy

#import re

#ff=subplot(111)


parser = OptionParser()
parser.add_option("-f", "--file",
    help="input file",
     metavar="FILE", default="all.csv", dest="file")

parser.add_option("-s", "--day",
    help="step to define a change in current",
    metavar="STEP", default="0.015",dest="step")


parser.add_option("-r", "--reqsteps",
        help="required steps  to become a plateau",
            metavar="RESTEPS", default="10",dest="reqsteps")



parser.add_option("-o", "--output",
                help="output file",
                metavar="out", default="out.dat", dest="out")

parser.add_option("-l", "--label",
                help="label for output files",
                metavar="LABEL", default="", dest="label")




def FindPlateaus(t,d,step,reqsteps):
  plts=[]
  AveRms=[]
  p0=t[0]
  p1="Not yet found"
  d0=float(d[0][0])
  i0=0
  i1=0
  for i in range(1,len(d)):
    if p1 != "Not yet found":
      i0=i
      p0=t[i]
      p1="Not yet found"
      d0=float(d[i][0])
    if (abs(float(d[i][0])-d0)>step or i==len(d)-1) and p1=="Not yet found":
      p1=t[i]
      if i-i0>reqsteps:  # This is a plateau
        plts.append([p0,p1])
        ave, rms, data = extractnew(t, d, p0,p1)
        
        AveRms.append([ave,rms])
      else:   # Things might be still changing this is not a plateau
        i0=i
        p0=t[i]
        p1="Not yet found"
        d0=float(d[i][0])

  return plts,AveRms



(options, args) = parser.parse_args()


data=parseout(options.file)

vars= data['datavars']
N=len(vars)
q1s= [string for string in vars if "RTQX1" in string]
Nq1=len(q1s)

print "Total number of vars: ", N
print "Vars for Q1 currents: ", Nq1

vB2x="LHC.BQBBQ.UA43.FFT1_B2:EIGEN_FREQ_1"
vB2y="LHC.BQBBQ.UA43.FFT1_B2:EIGEN_FREQ_2"
vB1x="LHC.BQBBQ.UA47.FFT1_B1:EIGEN_FREQ_1"
vB1y="LHC.BQBBQ.UA47.FFT1_B1:EIGEN_FREQ_2"




fout=open(options.out,"w")
print >>fout, "Day0  Hour0  Day1  Hour1 Current CurrentRms B2Qx  err  B2Qy  err   B1Qx   err   B1Qy   err "

refcu={}
refcu["RPMBC.UJ14.RTQX1.L1:I_MEAS"]=2.16190323333  
refcu["RPMBC.UJ16.RTQX1.R1:I_MEAS"]=-1.45739894263 
refcu["RPMBC.UJ56.RTQX1.R5:I_MEAS"]=-0.974910388412
refcu["RPMBC.USC55.RTQX1.L5:I_MEAS"]=-1.143901116

dd={}
data2show=[]
pltsDict={}
for iq in q1s:
  pltsDict[iq]=[]
  t=data[iq][0]
  d=data[iq][1]   # attention d elements are  lists

  plts, AveRms =FindPlateaus(t,d,float(options.step), int(options.reqsteps))
  print "Steps in ",iq,"    ->", len(plts)
  #print plts
  #print AveRms
  pltsDict[iq].append(plts)
  print >>fout, iq
  out2=open(iq+options.label,"w")
  print >>out2, "* Day0  Hour0  Day1  Hour1 Current CurrentRms B2Qx  B2Qxerr  B2Qy  B2Qyerr   B1Qx   B1Qxerr   B1Qy   B1Qyerr "
  print >>out2, "$ %s  %s  %s  %s %le %le %le  %le  %le  %le   %le   %le   %le   %le "

    
  
  for i in range(len(plts)):
    p0=plts[i][0]
    p1=plts[i][1]
    B2x = extract(data[vB2x][0], data[vB2x][1], p0,p1)
    B2y = extract(data[vB2y][0], data[vB2y][1], p0,p1)
    B1x = extract(data[vB1x][0], data[vB1x][1], p0,p1)
    B1y = extract(data[vB1y][0], data[vB1y][1], p0,p1)
    
    if  ("1:" in iq or "5:" in iq):
      if abs(refcu[iq]-AveRms[i][0])>0.4:
        print "selected",iq, AveRms[i][0]
        dd['name']=iq
        dd['current']=AveRms[i][0]
        dd['time0']=plts[i][0]
        dd['time1']=plts[i][1]
        dd['data']=[B2x[2],B2y[2],B1x[2],B1y[2]] # Send copies!
        data2show.append(dd.copy())
    print >>fout, plts[i][0],plts[i][1], AveRms[i][0], AveRms[i][1], B2x[0], B2x[1], B2y[0], B2y[1], B1x[0], B1x[1], B1y[0], B1y[1]
    print >>out2, plts[i][0],plts[i][1], AveRms[i][0], AveRms[i][1], B2x[0], B2x[1], B2y[0], B2y[1], B1x[0], B1x[1], B1y[0], B1y[1]

  out2.close()
fout.close()



def clean(d,l0,l1):
  return [x for x in d if (x<l1 and x>l0)]

def key(event):
  global case, data2show, PLANES, planecase
  #print "key", event.key, "figcanvas",event.canvas.figure
  if event.key=="c": close()
  if event.key=="n":
    if planecase==3:
      planecase=0
      case=case+1
    else:
      planecase=planecase+1
    clf()
    figure(1).canvas.set_window_title(PLANES[planecase])
    print "case:",case, "planecase",planecase 
    #data2show[case]['data'][planecase][0]
    hist(data2show[case]['data'][planecase],100)


def click (event):
    global lim0,lim1, data2show, fout, PLANES, planecase,case
    #print "you clicked" , event.x, event.y, event.xdata, event.ydata, event.button
    plot([event.xdata, event.xdata],[0, event.ydata])
    if event.button==1:
      lim0=event.xdata
      print "lim0 assigned" , event.xdata
    if event.button==3:
      lim1=event.xdata
      print "lim1 assigned" , event.xdata,"lim0 is", lim0
      newdata=clean(data2show[case]['data'][planecase], lim0, lim1)
      print average(newdata), std(newdata), PLANES[planecase]
      print >>fout, data2show[case]['name'],data2show[case]['time0'], data2show[case]['current'],   average(newdata), std(newdata), PLANES[planecase]
      
      event.key="n"
      key(event)





fout=open("cleanTunes.dat","w")
PLANES=["B2H","B2V","B1H","B1V"]
planecase=0
case=0
figure(1).canvas.set_window_title(PLANES[planecase])
#print data2show[case]['data'][planecase]
hist(data2show[case]['data'][planecase],100)
connect ( "button_press_event" , click )
connect ( "key_press_event" , key )

show()
fout.close()





