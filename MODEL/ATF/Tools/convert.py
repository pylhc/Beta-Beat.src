#
##
### Converter for ATF by Glenn Vanbacinckhove
##
#

from metaclass import twiss
from string import *
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", 
		 help="Filename that has to be converted",
		 metavar="FILE", default="null",dest="FILE")
parser.add_option("-m", "--model", 
		 help="Which bunch to use",
		 metavar="MODEL", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/ATF/",dest="MODEL")
parser.add_option("-o", "--output", 
		 help="Output path",
		 metavar="PATH", default="./",dest="PATH")

(options, args) = parser.parse_args()


# model
model=twiss(options.MODEL+"twiss.dat")

# tbt-data
tbt=open(options.FILE,'r')

# writing
filefile=open(options.PATH+"/"+os.path.basename(options.FILE)+".sdds.new","w")

lines=tbt.readlines()

monitors=[]
tbtdata=[]

# loading tbt into memory
for line in lines:

    if "n" in line: # adding header info

        print >> filefile,"#",line.rstrip('\n')


    elif "b" in line: # name of the monitors
        monitorline=line.split(" ")
        for monitor in monitorline:
            monitors.append(monitor.rstrip('\n'))
            
    else: # adding tbt

        tbtdata.append(line.rstrip('\n'))
            
# names of model
names=model.NAME
count=0
for mon in monitors:
    
    skip=0
    
    try:
       index=model.indx['M.'+mon[2:4]]
    except:
        print "WARNING : Monitor ",'M.'+mon[2:4]," not found will NOT print the tbt for this monitor"
        skip=1

    if skip==0:
        pos=model.S[index]
        plane=mon[4]

        if plane=="h":
            x=""
            for tbtx in tbtdata:
                xx=tbtx.split(" ")
                x=x+" "+xx[count]
                
            print >> filefile,"0",mon.replace("h",""),pos,x
        else:
            y=""
            for tbty in tbtdata:
                yy=tbty.split(" ")
                y=y+" "+yy[count]
            print >> filefile,"1",mon.replace("v",""),pos,y

    count=count+1


filefile.close()
