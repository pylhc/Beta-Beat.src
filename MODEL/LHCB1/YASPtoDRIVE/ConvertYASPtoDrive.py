from YASPmetaclassLHC import *
from optparse import OptionParser
import sys, gzip, os
from metaclass import *

######## optionparser
parser = OptionParser()
parser.add_option("-d", "--dir",
                help="which directory",
                metavar="DIR", default="./",dest="dir")
parser.add_option("-m", "--model",
                help="model location",
                metavar="TWISS", default="./",dest="twiss")
parser.add_option("-f", "--file",
                help="file to unzip",
                metavar="FILE", default="./",dest="file")

(options, args) = parser.parse_args()

###### main


dire=options.dir
output1=dire+"/multiturn1.dat"
output2=dire+"/multiturn2.dat"
fout1=open(output1, "w")
fout2=open(output2, "w")


print "Unzipping the files"
os.system("unzip "+options.file)

files = os.listdir(dire)
files = filter(lambda x: 'Dataset' in x    , files)
turns = len(files)

print turns, "files in dir:", dire

turnFiles={}

ListFiles=[]
for el in files:
    turn=int(el.split('-')[1])
    #print el, turn
    turnFiles[turn]=YASPtwiss(el)
    ListFiles.append(turnFiles[turn])
    


lhcb1=twiss(options.twiss+'/twiss.dat')
lhcb2=twiss(options.twiss+'/twiss.dat')



fout1.write('title\n')
fout2.write('title\n')

def intersectYASP(ListOfFile): 
	'''Pure intersection of all bpm names in all files '''
	if len(ListOfFile)==0:
		print "Nothing to intersect!!!!"
		sys.exit()
	h=ListOfFile[0].HNAME
	for b in ListOfFile:
		h=filter(lambda x: x in h   , b.HNAME)
        v=ListOfFile[0].VNAME
	for b in ListOfFile:
		v=filter(lambda x: x in v   , b.VNAME)
        
	
	return h, v

    

h , v = intersectYASP(ListFiles)

for el in h:
    if '.B1' in el:
        lhc=lhcb1
        f=fout1
    elif '.B2' in el:
        lhc=lhcb2
        f=fout2
    else:
        print "BPM ", el, "Not B1 Not B2 !!!!!!!!!"
        sys.exit()
    s=lhc.S[lhc.indx[el]]
    #print el, s
    f.write("0 "+el+" "+str(s)+" ")
    for i in range(1,turns+1):
        tF=turnFiles[i]
        f.write(str(tF.HCO[tF.Hindx[el]]*1000)+" ") #The *1000 is to convert to mum
    f.write("\n")



for el in v:
    if '.B1' in el:
        lhc=lhcb1
        f=fout1
    elif '.B2' in el:
        lhc=lhcb2
        f=fout2
    else:
        print "BPM ", el, "Not B1 Not B2 !!!!!!!!!"
        sys.exit()
    s=lhc.S[lhc.indx[el]]
    #print el, s
    f.write("1 "+el+" "+str(s)+" ")
    for i in range(1,turns+1):
        tF=turnFiles[i]
        f.write(str(tF.VCO[tF.Vindx[el]]*1000)+" ") #The *1000 is to convert to mum
    f.write("\n")


fout1.close()

fout2.close()



print "Removing unzipped data sets"
for b in files:

    print "removing",dire+"/"+b
    os.system("rm "+dire+"/"+b)

    

print "Finished !!"

    

        
    









