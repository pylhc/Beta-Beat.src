
###### imports
from optparse import OptionParser
from metaclass import twiss	
import os,sys,commands
from math import sqrt,cos,sin,pi,atan2
from datetime import date
from linreg import *


##
# YIL changes v 3.1:
#  - Cleaned macro writer in madcreator
#  - modifiers.madx now taken from options.twiss folder
#    if not found in translator[PATH]
#

###### optionparser
parser = OptionParser()
# general
parser.add_option("-m", "--model",
                help="twiss file to use",
                metavar="twiss", default="./", dest="twiss")
parser.add_option("-o", "--output",
                help="output path, where to store the results",
                metavar="output", default="./", dest="output")
parser.add_option("-b", "--beta",
                help="where beta-beat is stored",
                metavar="brc", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/", dest="brc")
parser.add_option("-t", "--technic",
                help="Which technique to use",
                metavar="technic", default="./", dest="technic")
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC",
                metavar="accel", default="LHCB1",dest="accel")

(options, args) = parser.parse_args()





## ############
#functions
## ############


#####
def madcreator(inifile,madfile,dpps):

	linesini=open(inifile,"r").readlines()
	linesmad=open(madfile+"/job.twiss_chrom.madx.macro","r").read()
	
	translator={}
	for line in linesini:
		li=line.split("=")
		translator[li[0].split()[0]]=li[1].split()[0]

	# creating the DPP
	dppstring=''
	dppstring_ac=''
	for dpp in dpps:
		if (os.path.exists(translator['PATH']+'/twiss_'+str(dpp)+'.dat')==False):
			dppstring=dppstring+'twiss, chrom,sequence='+translator['ACCEL']+', deltap='+str(dpp)+', file="'+translator['PATH']+'/twiss_'+str(dpp)+'.dat";\n'
			dppstring_ac=dppstring_ac+'twiss, chrom,sequence='+translator['ACCEL']+', deltap='+str(dpp)+', file="'+translator['PATH']+'/twiss_'+str(dpp)+'_ac.dat";\n'

	translator['DPP']=dppstring
	translator['DP_AC_P']=dppstring_ac

	for testpath in [translator['PATH'],options.twiss]:
		_tmpmod=os.path.join(testpath,'modifiers.madx')
		if os.path.isfile(_tmpmod):
			print "INFO: Using",_tmpmod
			translator['MODIFIERS']=_tmpmod
			break
	if 'MODIFIERS' not in translator:
		raise ValueError("Cannot find modifiers.madx")

	if(dppstring!=''):
		print "Creating madx"
		filetoprint=open(translator['PATH']+"/job.chrom.madx","w")


		#changing variables
		filetoprint.write(linesmad % translator)

		filetoprint.close()
		print "Running madx"
		os.system('madx < '+translator['PATH']+'/job.chrom.madx')



	else:
		print "No need to run madx"
			
###running getllm
def append(files):
	filestring="empty"

	for filee in files:
		filestring=filestring+","+filee

	return filestring.replace("empty,","")
		
def rungetllm(twissfile,accel,technic,files,outputpath,bsrc,dpp):

	VERSION="/GetLLM/GetLLM_V2.35.py"

	command="/usr/bin/python "+bsrc+VERSION+" -a "+accel+" -m "+twissfile+" -o "+outputpath+" -t "+technic+" -f "+append(files)

	print "Will run getllm for ",dpp, command
	
	os.system(command)
	print "GetLLM finished"

	os.system('cp '+outputpath+'/getbetax.out '+outputpath+'/getbetax_'+str(dpp)+'.out ')
	os.system('cp '+outputpath+'/getbetay.out '+outputpath+'/getbetay_'+str(dpp)+'.out ')
	os.system('cp '+outputpath+'/getampbetax.out '+outputpath+'/getampbetax_'+str(dpp)+'.out ')
	os.system('cp '+outputpath+'/getampbetay.out '+outputpath+'/getampbetay_'+str(dpp)+'.out ')
	os.system('cp '+outputpath+'/getcouple.out '+outputpath+'/getcouple_'+str(dpp)+'.out ')	
	os.system('cp '+outputpath+'/getbetax_free.out '+outputpath+'/getbetax_free_'+str(dpp)+'.out ')
	os.system('cp '+outputpath+'/getbetay_free.out '+outputpath+'/getbetay_free_'+str(dpp)+'.out ')
	os.system('cp '+outputpath+'/getcouple_free.out '+outputpath+'/getcouple_free_'+str(dpp)+'.out ')		

	

##### for chromatic
# model intersect
def modelIntersect(expbpms, model):

        bpmsin=[]
        for bpm in expbpms:
                try:
                        check=model.indx[bpm[1].upper()]
                        bpmsin.append(bpm)
                except:
                        print bpm, "Not in Model"
        if len(bpmsin)==0:
                print "Zero intersection of Exp and Model"
                print "Please, provide a good Dictionary"
                print "Now we better leave!"
                sys.exit()
        return bpmsin


#intersect
def intersect(ListOfFile): 
        '''Pure intersection of all bpm names in all files '''
        if len(ListOfFile)==0:
                print "Nothing to intersect!!!!"
                sys.exit()
        z=ListOfFile[0].NAME
        for b in ListOfFile:
                z=filter(lambda x: x in z   , b.NAME)
        #SORT by S
        result=[]
        x0=ListOfFile[0]
        for bpm in z:
                result.append((x0.S[x0.indx[bpm]], bpm))

        result.sort()
        return result

#linreg
def dolinregbet(filetoprint,listx,listy,bpms,plane,value,zero,twiss):
    for bpm in bpms:
	    el=bpm[1]
	    sloc=bpm[0]
	    indx=[]
	    b=[]
	    a=[]
	    bm=[]
	    am=[]
	    if "H" in plane:
		    beta0=zero.BETX[zero.indx[el]]
		    alfa0=zero.ALFX[zero.indx[el]]
		    alfa0err=zero.STDALFX[zero.indx[el]]
		    
		    beta0m=twiss.BETX[twiss.indx[el]] 
		    alfa0m=twiss.ALFX[twiss.indx[el]]
		    
		    wmo=twiss.WX[twiss.indx[el]]
		    pmo=twiss.PHIX[twiss.indx[el]]
	    else:
		    
		    beta0=zero.BETY[zero.indx[el]]
		    alfa0=zero.ALFY[zero.indx[el]]
		    alfa0err=zero.STDALFY[zero.indx[el]]

		    beta0m=twiss.BETY[twiss.indx[el]]
		    alfa0m=twiss.ALFY[twiss.indx[el]]
		    
		    wmo=twiss.WY[twiss.indx[el]]
		    pmo=twiss.PHIY[twiss.indx[el]]		    
	    for dpp in listx:
		    file=listy[dpp]
		    ix=file.indx[el]
		    indx.append(ix)
		    if "H" in plane:	
			    b.append(file.BETX[ix])
			    a.append(file.ALFX[ix])

			    bm.append(file.BETXMDL[file.indx[el]])
			    am.append(file.ALFXMDL[file.indx[el]])			    
		    else:
			    b.append(file.BETY[ix])
			    a.append(file.ALFY[ix])

			    bm.append(file.BETYMDL[file.indx[el]])
			    am.append(file.ALFYMDL[file.indx[el]])				    

	    bfit=linreg(listx, b)
	    afit=linreg(listx, a)

	    bfitm=linreg(listx, bm)
	    afitm=linreg(listx, am)

	    # measurement
	    dbb=bfit[0]/beta0
	    dbberr=bfit[3]/beta0
	    da=afit[0]
	    daerr=afit[3]
	    A=dbb
	    Aerr=dbberr
	    B=da-alfa0*dbb
	    Berr=sqrt(daerr**2 + (alfa0err*dbb)**2 + (alfa0*dbberr)**2)
	    w=sqrt(A**2+B**2)
	    werr=sqrt( (Aerr*A/w)**2 + (Berr*B/w)**2  )
	    phi=atan2(B,A)/2./pi
	    phierr=1./(1.+(A/B)**2)*sqrt( (Aerr/B)**2 + (A/B**2*Berr)**2)/2./pi

	    #model
	    dbbm=bfitm[0]/beta0m
	    dbberrm=bfitm[3]/beta0m
	    dam=afitm[0]
	    daerrm=afitm[3]
	    Am=dbbm
	    Aerrm=dbberrm
	    Bm=dam-alfa0m*dbbm
	    Berrm=sqrt(daerrm**2 + (alfa0m*dbberrm)**2)
	    wm=sqrt(Am**2+Bm**2)
	    werrm=sqrt( (Aerrm*Am/wm)**2 + (Berrm*Bm/wm)**2  )
	    phim=atan2(Bm,Am)/2./pi
	    phierrm=1./(1.+(Am/Bm)**2)*sqrt( (Aerrm/Bm)**2 + (Am/Bm**2*Berrm)**2)/2./pi
	    
	    
	    print >>filetoprint, el, sloc,  dbb, dbberr, da, daerr, w, werr, wmo,phi, phierr,pmo, dbbm,dbberrm,dam,daerrm,wm, werrm,phim,phierrm
    filetoprint.close()

### for coupling
# get det(C)
def getC(couplefile,name):

	
	f1001R=couplefile.F1001R[couplefile.indx[name]]
	f1001I=couplefile.F1001I[couplefile.indx[name]]	
	f1010R=couplefile.F1010R[couplefile.indx[name]]
	f1010I=couplefile.F1010I[couplefile.indx[name]]	

	down=4*((complex(f1001R,f1001I))-(complex(f1010R,f1010I)))
        c=1-(1/(1+down))

	cr=c .real
	ci=c .imag

	return cr,ci

# linreg for coupling
def dolinregCoupling(couplelist,bpms,dpplist,filetoprint,model):


	for bpm in bpms:

		name=bpm[1]
		s=bpm[0]
		
		a=[]
		br=[]
		bi=[]
	
		for dpp in dpplist:

			cr,ci=getC(couplelist[dpp],name)

		        a.append(dpp)
		        br.append(cr)
		        bi.append(ci)			
		
			
		else:

			fitr=linreg(a,br)
			fiti=linreg(a,bi)
			
			c=abs(complex(fitr[0],fiti[0]))
			e=abs(complex(fitr[3],fiti[3]))

			print >> filetoprint,name,s,c,e,"0"


		
			
	

	
	

## ##############
#main
## ##############

files=options.files.split(",")
outputpath=options.output
bsrc=options.brc
accel=options.accel
technic=options.technic


dpplist=[]
fileslist={}

for file in files:

	datax=twiss(file+"_linx")
	datay=twiss(file+"_liny")
	dppx=datax.DPP
	dppy=datay.DPP

	if dppx!=dppy:
		print "Discrepancy between horizontal and vertical => ",dppx,dppy
		print "System exit"
		sys.exit()
	else:
		dpp=dppx/1.0

#	if abs(dpp)==0.0004:
	#	print "ignoring"

	#else:
	
	if dpp not in dpplist:
		print "Adding dpp",dpp
		dpplist.append(dpp)
		fileslist[dpp]=[file]
	else:
		templist=fileslist[dpp]
		templist.append(file)
		print "The length of the list is ",len(templist)," for DPP ",dpp
		fileslist[dpp]=templist

if 0 not in dpplist:
	print "NO DPP=0.0"
	sys.exit()

madcreator(outputpath+"/super.ini",options.brc+"/MODEL/LHCB/model/",dpplist)
print "All models are created"
for dpp in dpplist:
	files=fileslist[dpp]
	rungetllm(outputpath+"/twiss_"+str(dpp)+".dat",accel,technic,files,outputpath,bsrc,dpp)
	#rungetllm(outputpath+"/twiss_0.0.dat",accel,technic,files,outputpath,bsrc,dpp)


##adding data
betalistx={}
betalisty={}
couplelist={}
betalistxf={}
betalistyf={}
couplelistf={}

listx=[]
listxf=[]
listy=[]
listyf=[]
listc=[]
listcf=[]

try:
	twiss(outputpath+'/getbetax_free_'+str(dpp)+'.out')
	freeswitch=1
except:
	freeswitch=0	


for dpp in dpplist:




        print "Loading driven data for ",dpp
        betx=twiss(outputpath+'/getbetax_'+str(dpp)+'.out')
        bety=twiss(outputpath+'/getbetay_'+str(dpp)+'.out')
        couple=twiss(outputpath+'/getcouple_'+str(dpp)+'.out')
	#couple=twiss(outputpath+'/getbetay_'+str(dpp)+'.out')
        betalistx[dpp]=betx
        betalisty[dpp]=bety
        couplelist[dpp]=couple

	if float(dpp)==0.0:
		zerobx=betx
		zeroby=bety

	listx.append(betx)
	listy.append(bety)
	listc.append(couple)
	modeld=twiss(options.twiss+"/twiss.dat")

        #try:
	if freeswitch==1:
                print "Loading free data"
                freeswitch=1
		print 'getbetax_free_'+str(dpp)+'.out'
                betxf=twiss(outputpath+'/getbetax_free_'+str(dpp)+'.out')
                betyf=twiss(outputpath+'/getbetay_free_'+str(dpp)+'.out')
                couplef=twiss(outputpath+'/getcouple_free_'+str(dpp)+'.out')
                betalistxf[dpp]=betxf
                betalistyf[dpp]=betyf
                couplelistf[dpp]=couplef
		listxf.append(betxf)
		listyf.append(betyf)
		listcf.append(couplef)
		modeld=twiss(options.twiss+"/twiss_ac.dat")
		modelf=twiss(options.twiss+"/twiss.dat")
		if float(dpp)==0.0:
			zerobxf=betalistxf[dpp]
			zerobyf=betalistyf[dpp]	
			

        #except:
         #       print "No free data"

#
# driven beta
#

print "Driven beta"

#H
filefile=open(outputpath+"/chrombetax.out","w")
print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WX","WXERR","WMO","PHIX", "PHIXERR","PHIM", "dbbR", "dbberrR", "dalfaR", "daerr","WXR","WXERRR","PHIXR", "PHIXERRR"
print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le %le  %le %le  %le %le  %le %le %le  %le %le  %le"

bpms=intersect(listx)
bpms=modelIntersect(bpms,modeld)
dolinregbet(filefile,dpplist,betalistx,bpms,"H","beta",zerobx,modeld)
filefile.close()

#V
filefile=open(outputpath+"/chrombetay.out","w")
print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WY", "WYERR","WYM","PHIY",  "dbbR", "dbberrR", "dalfaR", "daerr","PHIYERR","PHIM", "WYR","WYERRR","PHIYR", "PHIYERRR"
print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le %le %le %le  %le %le %le %le  %le %le %le"

bpms=intersect(listy)
bpms=modelIntersect(bpms,modeld)
dolinregbet(filefile,dpplist,betalisty,bpms,"V","beta",zeroby,modeld)
filefile.close()

print "Driven beta finished"

#
# driven coupling
#
print "Driven coupling"

filefile=open(outputpath+"/chromcoupling.out","w")
print >>filefile,"NAME S CHROMCOUPLE  CHROMe  CHROMMDL"
print >>filefile,"%s   %le  %le       %le     %le"

bpms=intersect(listc)
bpms=modelIntersect(bpms,modeld)

dolinregCoupling(couplelist,bpms,dpplist,filefile,modeld)
filefile.close()


print "Driven coupling finished"
filefile.close()

if freeswitch==1:
  #
  # free beta
  #
  print "Free beta"
  #H
  filefile=open(outputpath+"/chrombetax_free.out","w")
  print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WX","WXERR","WMO","PHIX", "PHIXERR","PHIM",  "dbbR", "dbberrR", "dalfaR", "daerr","WXR","WXERRR","PHIXR", "PHIXERRR"
  print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le %le %le %le  %le %le %le %le  %le"

  bpms=intersect(listxf)
  bpms=modelIntersect(bpms,modelf)
  dolinregbet(filefile,dpplist,betalistxf,bpms,"H","beta",zerobxf,modelf)
  filefile.close()

  #V
  filefile=open(outputpath+"/chrombetay_free.out","w")
  print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WY", "WYERR","WYM","PHIY", "PHIYERR","PHIM",  "dbbR", "dbberrR", "dalfaR", "daerr","WYR","WYERRR","PHIYR", "PHIYERRR"
  print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le %le %le %le  %le %le %le %le  %le"

  bpms=intersect(listyf)
  bpms=modelIntersect(bpms,modelf)
  dolinregbet(filefile,dpplist,betalistyf,bpms,"V","beta",zerobyf,modelf)
  filefile.close()

  print "Free beta finished"

  #
  # free coupling
  #
  print "Free coupling"

  filefile=open(outputpath+"/chromcoupling_free.out","w")
  print >>filefile,"NAME S CHROMCOUPLE  CHROMe  CHROMMDL"
  print >>filefile,"%s   %le  %le       %le     %le"

  bpms=intersect(listcf)
  bpms=modelIntersect(bpms,modelf)

  dolinregCoupling(couplelistf,bpms,dpplist,filefile,modelf)
  filefile.close()


  print "Free coupling finished"
  filefile.close()


sys.exit()

