import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD/")
from pylab import *
from Timber import *
from optparse import OptionParser
import numpy, sys,  metaclass


#import re

#ff=subplot(111)


parser = OptionParser()
parser.add_option("-f", "--file",
    help="input file with coupling and tunes",
     metavar="FILE", default="testb2ir1.csv", dest="file")



parser.add_option("-p", "--plateaus",
        help="File containing timestamps for plateaus",
            metavar="PFILE", default="testb2ir1plat.tfs",dest="pfile")



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

plateaus=metaclass.twiss(options.pfile)

vars= data['datavars']
N=len(vars)

vB2x="LHC.BQBBQ.UA43.FFT1_B2:EIGEN_FREQ_1"
vB2y="LHC.BQBBQ.UA43.FFT1_B2:EIGEN_FREQ_2"
vB2c="LHC.BQBBQ.UA43.FFT1_B2:COUPLING_ABS"
#vB1x="LHC.BQBBQ.UA47.FFT1_B1:EIGEN_FREQ_1"
#vB1y="LHC.BQBBQ.UA47.FFT1_B1:EIGEN_FREQ_2"



#print vars
print "Total number of vars: ", N
print "Number of plateaus: ", len(plateaus.Date)


# Defining plts
plts=[]
AveRms=[]
for i in range(1,len(plateaus.Date)):
  plts.append([plateaus.Date[i-1]+" "+plateaus.Time[i-1], plateaus.Date[i]+" "+plateaus.Time[i]] )
  #print plts[-1]
  AveRms.append([plateaus.Value[i-1], 0])


fout=open(options.out,"w")



print "Steps     ->", len(plts)
  #print plts
  #print AveRms
  


print >>fout, "* Day0  Hour0  Day1  Hour1 Current CurrentRms B2Qx  B2Qxerr  B2Qy  B2Qyerr B2c-  B2c-err"
print >>fout, "$ %s  %s  %s  %s %le %le %le  %le  %le  %le   %le   %le   %le   %le "

    
dd={}
data2show=[]
pltsDict={}

 
for i in range(len(plts)):
    p0=plts[i][0]
    p1=plts[i][1]
    B2x = extractnew(data[vB2x][0], data[vB2x][1], p0,p1)
    B2y = extractnew(data[vB2y][0], data[vB2y][1], p0,p1)
    B2c = extractnew(data[vB2c][0], data[vB2c][1], p0,p1)
    #B1x = extractnew(data[vB1x][0], data[vB1x][1], p0,p1)
    #B1y = extractnew(data[vB1y][0], data[vB1y][1], p0,p1)
    
    
      
        
    dd['name']=options.file
    dd['current']=AveRms[i][0]
    dd['time0']=plts[i][0]
    dd['time1']=plts[i][1]
    dd['data']=[B2x[2],B2y[2],B2c[2]] # Send copies!
    data2show.append(dd.copy())
    print i
    print >>fout, plts[i][0],plts[i][1], AveRms[i][0], AveRms[i][1], B2x[0], B2x[1], B2y[0], B2y[1], B2c[0], B2c[1]
    #print >>out2, plts[i][0],plts[i][1], AveRms[i][0], AveRms[i][1], B2x[0], B2x[1], B2y[0], B2y[1], B1x[0], B1x[1], B1y[0], B1y[1]

  
fout.close()



def clean(d,l0,l1):
  return [x for x in d if (x<l1 and x>l0)]

def limitset():
  global limset, setlim
  if setlim==1:
    print "limits set"
    limset=1
  else:
    print "no limit set"
    limset=0
  setlim=0
  return limset

def nextplot():
  global case, data2show, PLANES, planecase
  if planecase==2:
    planecase=0
    case=case+1
  else:
    planecase=planecase+1
  clf()
  figure(1).canvas.set_window_title(PLANES[planecase]+" #:"+str(case))
  print "case:",case, "planecase",planecase 
  #data2show[case]['data'][planecase][0]
  trimq,q,qer = loadtxt('dispQ.dat', delimiter=' ', unpack=True, usecols=[0,1,2])
  if case>0:
    trimc,c,cer = loadtxt('dispC.dat', delimiter=' ', unpack=True, usecols=[0,1,2])
  subplot(212)
  hist(data2show[case]['data'][planecase],100)
  subplot(222)
  plot1=errorbar(meandata[:,0],meandata[:,2],yerr=meandata[:,3], fmt='r.')
  plot2=errorbar(meandata[:,0],meandata[:,4],yerr=meandata[:,5], fmt='c.')
  plot3=errorbar(trimq, q, yerr=qer, fmt='k.')
  subplot(221)
  plot4=errorbar(meandata[:,0],meandata[:,6],yerr=meandata[:,7], fmt='g.')
  if case>0:
    plot5=errorbar(trimc, c, yerr=cer, fmt='k.')

def cleanplot():
  if planecase==0 or planecase==1:
    subplot(222)
    plot3=errorbar(data2show[case]['current'],average(newdatatemp),std(newdatatemp), fmt='k.')
  elif planecase==2:
    subplot(221)
    plot5=errorbar(data2show[case]['current'],average(newdatatemp),std(newdatatemp), fmt='k.')

def key(event):
  global case, data2show, PLANES, planecase, newdata
  if event.key=="c":
    close()
    sys.exit("cleaning was stoped by pressing c")
  if event.key=="x":
    print "plateau excluded"
    nextplot()
  elif event.key=="n":
    limitset()
    if limset==0:
      newdata=data2show[case]['data'][planecase]
      print "no cleaning performed"
    elif limset==1:
      newdata=clean(data2show[case]['data'][planecase], lim0, lim1)
    print average(newdata), std(newdata), PLANES[planecase]
    print >>fout, data2show[case]['name'],data2show[case]['time0'], data2show[case]['current'],   average(newdata), std(newdata), PLANES[planecase]
    dispQ=open("dispQ.dat","a")
    dispC=open("dispC.dat","a")
    if planecase==0 or planecase==1:
      print >>dispQ, data2show[case]['current'], average(newdata), std(newdata), PLANES[planecase]
    elif planecase==2:
      print >>dispC, data2show[case]['current'], average(newdata), std(newdata), PLANES[planecase]
    dispQ.close()
    dispC.close()
    
    nextplot()


def click (event):
    global setlim, lim0, lim1, data2show, fout, PLANES, planecase, case, newdatatemp
    setlim=0
    #print "you clicked" , event.x, event.y, event.xdata, event.ydata, event.button
    subplot(212).plot([event.xdata, event.xdata],[0, event.ydata])
    if event.button==1:
      lim0=event.xdata
      print "lim0 assigned" , event.xdata
    if event.button==3:
      lim1=event.xdata
      print "lim1 assigned" , event.xdata,"lim0 is", lim0
      setlim=1
      newdatatemp=clean(data2show[case]['data'][planecase], lim0, lim1)
      cleanplot()





#fname=open('./out.dat','r')
#fname=cbook.get_sample_data('msft.csv', asfileobj=False)
fout=open("cleanTunes.dat","w")
dispQ=open("dispQ.dat","w")
dispQ.close()
dispC=open("dispC.dat","w")
dispC.close()
setlim=0
PLANES=["B2H","B2V","B2C"]
planecase=0
case=0
meandata=loadtxt('out.dat',skiprows=1, comments='$', delimiter=' ', usecols=[4,5,6,7,8,9,10,11])
figure(1).canvas.set_window_title(PLANES[planecase])
ion()
print "Interactive mode of pylabe?", isinteractive()
subplot(212)
ploth=hist(data2show[case]['data'][planecase],100)
subplot(222)
plot1=errorbar(meandata[:,0],meandata[:,2],yerr=meandata[:,3], fmt='r.')
plot2=errorbar(meandata[:,0],meandata[:,4],yerr=meandata[:,5], fmt='c.')
subplot(221)
plot4=errorbar(meandata[:,0],meandata[:,6],yerr=meandata[:,7], fmt='g.')


#print data2show[case]['data'][planecase]
#figure(2).canvas.set_window_title("out.dat")
#plotfile(fname, ('Current','B2Qx','B2Qy'))
#plotfile(fname, (0,5,6))
connect ( "button_press_event" , click )
connect ( "key_press_event" , key )

show()
fout.close()





