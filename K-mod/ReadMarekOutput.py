from metaclass25 import *
import sys, os
from optparse import OptionParser


parser = OptionParser()
#parser.add_option("-a", "--accel",
#                help="Which accelerator: LHCB1 LHCB2 ",
#                metavar="ACCEL", default="LHCB1",dest="ACCEL")

parser.add_option("-f", "--file",
                help="Input file containg DKs and DQs as output from Marek's code",
                metavar="FILE", default="", dest="file")


parser.add_option("-m", "--modifiers",
                help="modifiers MADX file to define beta*",
                metavar="", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/modifiers.madx", dest="modifiers")


(options, args) = parser.parse_args()



f=twiss(options.file, "r")


dic={"Q1R1":"MQXA.1R1",
     "Q1L1":"MQXA.1L1",
     "Q1R2":"MQXA.1R2",
     "Q1L2":"MQXA.1L2",
     "Q1R5":"MQXA.1R5",
     "Q1L5":"MQXA.1L5",
     "Q1R8":"MQXA.1R8",
     "Q1L8":"MQXA.1L8"}


print "Looking for quads in the input"

IRlabel={"IR1":0, "IR2":0 ,"IR5":0, "IR8":0  }


f.DEVICE=f.NAME

if dic["Q1R1"] in f.DEVICE and dic["Q1L1"] in f.DEVICE:
    print "IR1 quads found"
    IRlabel["IR1"]=1
if dic["Q1R2"] in f.DEVICE and dic["Q1L2"] in f.DEVICE:
    print "IR2 quads found"
    IRlabel["IR2"]=1
if dic["Q1R5"] in f.DEVICE and dic["Q1L5"] in f.DEVICE:
    print "IR5 quads found"
    IRlabel["IR5"]=1
if dic["Q1R8"] in f.DEVICE and dic["Q1L8"] in f.DEVICE:
    print "IR8 quads found"
    IRlabel["IR8"]=1


for IR in IRlabel.keys():
    if IRlabel[IR]==1:
        print "Processing ", IR
        IRn=IR[-1:]
        QR=dic["Q1R"+IRn]
        QL=dic["Q1L"+IRn]
        iQR=f.DEVICE.index(QR)
        iQL=f.DEVICE.index(QL)
        kQR=f.deltaK[iQR]
        kQL=f.deltaK[iQL]
        if (kQR != kQL):
            print "Sorry, I do not like different k's in", IR
            sys.exit()
        print "# Beam 1 H and V"
        print "python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/MatchQ.py -a LHCB1 -I"+IRn+" -pQ1"+" -d"+str(f.dQH_B1[iQR])+","+str(f.dQH_B1[iQL])+" -k"+str(kQR), "-m ", options.modifiers
        print "python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/MatchQ.py -a LHCB1 -I"+IRn+" -pQ2"+" -d"+str(f.dQV_B1[iQR])+","+str(f.dQV_B1[iQL])+" -k"+str(kQR), "-m ", options.modifiers
        print "# Beam 2 H and V"
        print "python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/MatchQ.py -a LHCB2 -I"+IRn+" -pQ1"+" -d"+str(f.dQH_B2[iQR])+","+str(f.dQH_B2[iQL])+" -k"+str(kQR), "-m ", options.modifiers
        print "python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/MatchQ.py -a LHCB2 -I"+IRn+" -pQ2"+" -d"+str(f.dQV_B2[iQR])+","+str(f.dQV_B2[iQL])+" -k"+str(kQR), "-m ", options.modifiers

        
# Command MatchQ example
#python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/MatchQ.py -a LHCB1 -I5 -pQ1 -d0.0009,0.00105 -k1e-5


