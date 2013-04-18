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
output1=dire+"/multiturn1.dat.sdds"
output2=dire+"/multiturn2.dat.sdds"
fout1=open(output1, "w")
fout2=open(output2, "w")


#print "Unzipping the files"
#os.system("unzip "+options.file +" -d "+options.dir)

filefile=open(options.file,"r")
lines=filefile.readlines()

#step=5
#count=0

for line in lines:

    if "<measureddatafile>" in line:
        print "Found line to path"
        line=line.split('>')[1]
        line=line.split('<')[0]
        os.system("cp "+line+" "+options.dir)
        base=os.path.basename(line)
        os.system("gunzip "+base )
            #os.system("rm "+base.split('.gz')[0])
        count=count+1


files = os.listdir(dire)
files = filter(lambda x: 'ORBIT_LHCRING' in x    , files)
turns = len(files)
#turns=120  ###### quick fix!!!!!

print turns, "files in dir:", dire


turnFiles={}

ListFiles=[]
counter=0
for el in files:
    #turn=int(el.split('-')[1])
#    print el, turn
    if counter<len(files)+1:
#        print el,turn
        turnFiles[counter]=YASPtwiss(dire+"/"+el)
        ListFiles.append(turnFiles[counter])
        counter=counter+1

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
        count=0
	for b in ListOfFile:
            if len(b.HNAME)==0:
                break
            else:
                h=filter(lambda x: x in h   , b.HNAME)
                count=count+1
            #h=filter(lambda x: x in h   , b.HNAME)
        turns=count  
        v=ListOfFile[0].VNAME
	for b in ListOfFile:
            if len(b.VNAME)==0:
                break
            else:
                v=filter(lambda x: x in v   , b.VNAME)
            v=filter(lambda x: x in v   , b.VNAME)

                
        print "Number of turns that can be used ",str(turns)
	
	return h, v,turns

    

h , v ,turns= intersectYASP(ListFiles)



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
    for i in range(1,len(files)):
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
    for i in range(1,len(files)):
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

    

        
    









