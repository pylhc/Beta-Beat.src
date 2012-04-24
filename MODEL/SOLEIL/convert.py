#################################
#
# @ Glenn Vanbavinckhove
# convert data to coupling data
#
#################################

from metaclass import *
from optparse import OptionParser
from string import *

# getting options
parser = OptionParser()

parser.add_option("-b", "--bpm",
                help="bpm coupling",
                metavar="BPM", default="./BPMcoupling.tfs",dest="Bpm")
parser.add_option("-f", "--file",
                help="file that has to be converted",
                metavar="FILE", default="./",dest="File")
(options, args) = parser.parse_args()
# loading file
convertF=options.Bpm
conFile=twiss(convertF)
f = open(options.File,'r');x=f.readlines();f.close()
fdes=open(options.File+'.sdds.new','w')

# main
for i in range(len(x)-1):
    
     s = split(x[i])
     
     
     if len(s) >= 2:
         bpm=s[1]
        

         d=s[0]
         if d[0]=='#':
            for j in range(0,len(s)):
                if j==len(s)-1:
                    fdes.write(str(s[j])+' \n')
                else:
                    fdes.write(str(s[j])+' ')

         else:

             hbpmgain=conFile.HBPMGain[conFile.indx[bpm]]
             hbpmcoupling=conFile.HBPMCoupling[conFile.indx[bpm]]
            

             vbpmgain=conFile.VBPMGain[conFile.indx[bpm]]
             vbpmcoupling=conFile.VBPMCoupling[conFile.indx[bpm]]

            
             
             sx= split(x[i])
             sy= split(x[i+1])
         
             data=[]
            
             if sx[0]=='0':
                  if sx[1]!=sy[1]:
                       print "Converter is not ready for this BPM data !! sorry..."
                       sys.exit()
                  if  sx[0]=='1':
                       print "plane seems to be switched, converter not ready for this !! sorry..."
                       sys.exit()
                       
                  fdes.write(sx[0]+' ')
                  fdes.write(sx[1]+' ')
                  fdes.write(sx[2]+' ')
         
                  for j in range(0,len(sx)):

                      if j >=3:
         
                          xx=hbpmgain*float(sx[j])+hbpmcoupling*float(sy[j])
                          if j==len(sx)-1:
                              fdes.write(str(xx)+' \n')
                          else:
                              fdes.write(str(xx)+' ')
   

                  fdes.write(sy[0]+' ')
                  fdes.write(sy[1]+' ')
                  fdes.write(sy[2]+' ')
                  for j in range(0,len(sy)):
                         
                      if j >=3:
                          
                         yy=vbpmgain*float(sy[j])+vbpmcoupling*float(sx[j])
                         if j==len(sy)-1:
                             fdes.write(str(yy)+' \n')
                         else:
                             fdes.write(str(yy)+' ')

fdes.close()
