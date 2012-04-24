import sys
from optparse import OptionParser
from string import *
import re
#from os import system
import os

os.system('export LD_LIBRARY_PATH; LD_LIBRARY_PATH=/usr/X11R6/lib:/lib:/usr/lib:/usr/local/lib:/usr/local/ActiveTcl/lib:/usr/local/share/sybase/OCS-12_5/lib:/ride/release/X86/lib')


def cbpm(s): # some model convention problems    
    s=replace(s,"rbpm.","");s=replace(s,"bpm.","")
    s=replace(s,"b-g","g");s=replace(s,"y-g","g")
    s=replace(s,"-bhx","_bx");s=replace(s,"-bvx","_bx")
    if re.match('.*-bh1$',s):s=replace(s,"-bh1","_b1")
    if re.match(".*-bv1$",s):s=replace(s,"-bv1","_b1")
    if re.match(".*-bv3$",s):s=replace(s,"-bv3","_b3")
    if re.match(".*-bh3$",s):s=replace(s,"-bh3","_b3")
    
    s=replace(s,"-bh4","_b4");s=replace(s,"-bv4","_b4")
    s=replace(s,"-bh7","_b7");s=replace(s,"-bv7","_b7")
    s=replace(s,"-bh8","_b8");s=replace(s,"-bv8","_b8")
    s=replace(s,"-bh3.1","_b3.1");s=replace(s,"-bv3.1","_b3.1")
    s=replace(s,"-bh3.2","_b3.2");s=replace(s,"-bv3.2","_b3.2")
    s=replace(s,"-bh7.1","_b7.1");s=replace(s,"-bv7.1","_b7.1")
    s=replace(s,"-b","_b")

    s=replace(s,'\n',' ')
    
    if re.match("^bo10_b3$",s): s=replace(s,"bo10_b3","bo10_b3.1")
    if re.match("^bi1_b3$",s): s=replace(s,"bi1_b3","bi1_b3.1")
    return s.upper()
    #return s


parser = OptionParser()
parser.add_option("-a", "--accel",
    help="Accelerator: LHCB1 LHCB2 SPS RHIC  SOLEIL",
    metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-f", "--files",
    help="Files from analysis, separated by comma",
    metavar="TwissFile", default="0", dest="files")
(opt, args) = parser.parse_args()

files=opt.files.split(",")

print 'number of files to convert '+str(len(files))


for filee in files:
    
    print 'converting binary to ascii'
    # convert binary to ascii
    print '/ride/release/X86/bin/RhicBpmAnalysis  -a  temp '+filee
    os.system('/ride/release/X86/bin/RhicBpmAnalysis  -a  temp '+filee )
    os.system('chmod 777 temp')
    os.system('chmod 777 '+filee )
    os.system('mv temp '+filee)
    
    # taking namins conventions into account

    print 'starting name convertions'
    f=open(filee,'r')
    
    print 'converting file '+str(filee)

    lines=f.readlines()

    d=open('./temp','w')
    
    for i in lines:
      
        dd=cbpm(i)
        print >> d, dd
 
    

    d.close()

    os.system('mv temp '+filee)
