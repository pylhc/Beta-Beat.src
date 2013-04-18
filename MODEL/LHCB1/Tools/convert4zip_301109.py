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
parser.add_option("-o","--optics",
		help="optics to use",
		metavar="OPT", default="nominal.opt",dest="opt")
parser.add_option("-b","--beam",
		help="beam to name file",
		metavar="BEAM", default="B1",dest="beam")

(options, args) = parser.parse_args()

###### main


dire=options.dir
if options.beam=="B1":
    output1=dire+"/"+os.path.basename(options.file)+".sdds.new"
else:    
    output1=dire+"/multiturn1.dat.sdds"
    
if options.beam=="B2":
    output2=dire+"/"+os.path.basename(options.file)+".sdds.new"
else:
    output2=dire+"/multiturn2.dat.sdds"
fout1=open(output1, "w")
fout2=open(output2, "w")


print "Unzipping the files"
os.system("unzip "+options.file +" -d "+options.dir)

files = os.listdir(dire)
files = filter(lambda x: 'Dataset' in x    , files)
#turns = len(files)
turns=120  ###### quick fix!!!!!

print turns, "files in dir:", dire


turnFiles={}

ListFiles=[]
for el in files:
    turn=int(el.split('-')[1])
#    print el, turn
    if turn<turns+1:
#        print el,turn
        turnFiles[turn]=YASPtwiss(dire+"/"+el)
        ListFiles.append(turnFiles[turn])

#sys.exit()


lhcb1=twiss(options.twiss+'/LHCB1/'+options.opt+'/twiss.dat')
lhcb2=twiss(options.twiss+'/LHCB2/'+options.opt+'/twiss.dat')



fout1.write('# title\n')
fout2.write('# title\n')

def intersectYASP(ListOfFile): 
	'''Pure intersection of all bpm names in all files '''
	if len(ListOfFile)==0:
		print "Nothing to intersect!!!!"
		sys.exit()
	h=ListOfFile[0].HNAME
        #print len(h)
        count =0
	for b in ListOfFile:
                #print b
                #sys.exit()
                print len(b.HNAME),count
		h=filter(lambda x: x in h   , b.HNAME)
                count=count+1
        #sys.exit()
        v=ListOfFile[0].VNAME
        #print len(h)
        #sys.exit()
	for b in ListOfFile:
		v=filter(lambda x: x in v   , b.VNAME)
        print "the length of v and h ",str(len(h)),str(len(v))
	
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

    

        
    









