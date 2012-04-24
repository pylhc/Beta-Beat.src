# This program transforms turn-by-turn data into TBT with correction factor for the non-linearity of the BPM (xobs=xr-k*xr**3)
#
#@ author : Rogelio Garcia Tomas, Masamitsu Aiba and Glenn Vanbavinckhove
#
#@ date : 02/04/09


#from roots import *
from metaclass import *
from rhicdata import *
import math, cmath,os
from optparse import *


# optionparser
parser = OptionParser()
parser.add_option("-f", "--file base",
                      help="Path to file with twiss file",
                      metavar="FILE", default="./",dest="FILE")

parser.add_option("-k", "--kick",
                      help="Path to file with twiss file",
                      metavar="KICK", default="6",dest="KICK")
parser.add_option("-l", "--lin_x",
                      help="Path to file with twiss file",
                      metavar="LIN", default="./",dest="LIN")
parser.add_option("-o", "--output",
                      help="place to put the result",
                      metavar="OUT", default="./",dest="OUT")
parser.add_option("-c", "--calibration_file",
                      help="calibration file",
                      metavar="CAL", default="./",dest="CAL")
(opt, args) = parser.parse_args()

#
#Command => python CubicReconstruction.new.py -f ../blabla.txt -k 250 -l ../blabla.txt_linx -o ../here/
# -c ./calibration.tfs


#V=opt.KICK
lin=opt.LIN
#lx=twiss(lin)
filee=opt.FILE
cali=twiss(opt.CAL)
data=rhicdata(filee)

#xmm = 1.7  # x calibration factor
#zmm = 0.58 # z (y CERN) calibration factor
#bxref=13.83834969 # beta ref
Klin=average(array(cali.KLINX))
Knonlin=average(array(cali.KNONLINX))

print 'linear factor is'+str(Klin)
print 'non-linear factor is'+str(Knonlin)

#m=twiss('twiss.dat')

name=os.path.basename(filee)
f=open(opt.OUT+'/'+name,"w")
f.write(data.title)



for i in range(len(data.H)):
    name=data.H[i].name
    s=data.H[i].location

    #print name, V,  a_30, kapprox
    f.write("0 "+name+" "+str(s)+" ")
    for j in range(len(data.H[i].data)):
        xobs=data.H[i].data[j]


        xreal=(1/Klin[i])*xobs-(Knonlin[i]//Klin[i]**4)*xobs**3
        #xreal=Klin[i]*xobs-Knonlin[i]*xobs**3
  
        f.write(str(xreal)+" ")
    f.write("\n")

    
    f.write("1 "+name+" "+str(s)+" ")
    for j in range(len(data.V[i].data)):
        f.write(str(data.V[i].data[j])+" ")
    f.write("\n")

f.close()
