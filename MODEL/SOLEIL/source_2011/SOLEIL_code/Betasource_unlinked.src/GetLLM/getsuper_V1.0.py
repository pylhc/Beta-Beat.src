###### imports
from optparse import OptionParser
try:
	from metaclass import twiss
except:
	from metaclass25 import twiss	
import os,sys
from math import sqrt,cos,sin,pi,atan2
from datetime import date
from linreg import *




###### optionparser
parser = OptionParser()
parser.add_option("-f", "--files",
                help="Files to use (seperated by ,)",
                metavar="files", default="./", dest="files")
parser.add_option("-t", "--twiss",
                help="twiss file to use",
                metavar="twiss", default="./", dest="twiss")
parser.add_option("-o", "--output",
                help="output path, where to store the results",
                metavar="output", default="./", dest="output")
parser.add_option("-d", "--dpp",
                help="dpps to take",
                metavar="dpps", default="./", dest="dpps")

(options, args) = parser.parse_args()





## ############
#functions
## ############

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
	    if "H" in plane:
		    beta0=zero.BETX[zero.indx[el]]
		    alfa0=zero.ALFX[zero.indx[el]]
		    alfa0err=zero.STDALFX[zero.indx[el]]
		    wmo=twiss.WX[twiss.indx[el]]
		    pmo=twiss.PHIX[twiss.indx[el]]
	    else:
		    
		    beta0=zero.BETY[zero.indx[el]]
		    alfa0=zero.ALFY[zero.indx[el]]
		    alfa0err=zero.STDALFY[zero.indx[el]]
		    wmo=twiss.WY[twiss.indx[el]]
		    pmo=twiss.PHIY[twiss.indx[el]]		    
	    for dpp in listx:
		    file=listy[dpp]
		    ix=file.indx[el]
		    indx.append(ix)
		    if "H" in plane:	
			    b.append(file.BETX[ix])
			    a.append(file.ALFX[ix])
		    else:
			    b.append(file.BETY[ix])
			    a.append(file.ALFY[ix])			

	    bfit=linreg(listx, b)
	    afit=linreg(listx, a)
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
	    
	    print >>filetoprint, el, sloc,  dbb, dbberr, da, daerr, w, werr, wmo,phi, phierr,pmo

    filetoprint.close()

### for coupling
# get det(C)
def getC(couplefile,name):

	
	f1001=couplefile.F1001W[couplefile.indx[name]]
	f1010=couplefile.F1010W[couplefile.indx[name]]

	check=(1/4)+f1010**2

	if check==f1001**2:
		print "Skipping"
		skip="Yes"
		C=0
	else:
		part1= 4*(f1001**2-f1010**2)
		C=1-(1/(1+part1))
		skip="No"

	return C,skip

# linreg for coupling
def dolinregCoupling(couplelist,bpms,dpplist,filetoprint,model):


	for bpm in bpms:

		name=bpm[1]
		s=bpm[0]
		
		a=[]
		b=[]
	
		for dpp in dpplist:

			c,skip=getC(couplelist[dpp],name)

			if skip=="No":

				a.append(dpp)
				b.append(c)
			else:
				a.append(dpp)
				b.append(skip)
				

		if "Yes" in b:

			print "Skipping chromatic coupling calculation for ",bpm

		else:

			fit=linreg(a, b)

			print >> filetoprint,name,s,fit[0],fit[3],"0"


		
			
	

	
	

## ##############
#main
## ##############

files=options.files.split(",")
dpps=options.dpps.split(",")
model=twiss(options.twiss)
output=options.output

if len(dpps)!=len(files):
	print "Unequal input"
	sys.exit()

betalistx={}
listx=[]
betalisty={}
listy=[]
couplelist={}
couplel=[]


betalistxf={}
betalistyf={}
couplelistf={}

dpplist=[]


for count in range(len(files)):

	filee=files[count]
	dpp=float(dpps[count])

	print "Loading for ",dpp


	print "Loading driven data"
	betx=twiss(filee+'/getbetax.out')
	bety=twiss(filee+'/getbetay.out')
	couple=twiss(filee+'/getcouple.out')
	dpplist.append(float(dpp))
	
	listx.append(betx)
	betalistx[dpp]=betx
	listy.append(bety)
	betalisty[dpp]=bety
	couplel.append(couple)
	couplelist[dpp]=couple


	try:
		print "Loading free data"
		freeswitch=1
		betxf=twiss(filee+'/getbetax_free.out')
		betyf=twiss(filee+'/getbetay_free.out')
		couplef=twiss(filee+'/getcouple_free.out')
		betalistxf[dpp]=betxf
		betalistyf[dpp]=betyf
		couplelistf[dpp]=couplef
		
	except:
		print "No free data"


	

### finding dpp index
for i in range(len(dpps)):
	
    if float(dpps[i])==0.0:
	dpp=float(dpps[i])
        zeroi=i
        zerobx=betalistx[dpp]
        zeroby=betalisty[dpp]
        zerobxf=betalistxf[dpp]
        zerobyf=betalistyf[dpp]	
        print "Found dpp=0 for case:", files[i]
if zeroi<0:
    print "dpp=0 not found, exit"
    sys.exit()



#
# driven beta
#

print "Driven beta"

#H
filefile=open(output+"/chrombetax.out","w")
print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WX","WXERR","WMO","PHIX", "PHIXERR","PHIM"
print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le %le  %le"

bpms=intersect(listx)
bpms=modelIntersect(bpms,model)
dolinregbet(filefile,dpplist,betalistx,bpms,"H","beta",zerobx,model)
filefile.close()

#V
filefile=open(output+"/chrombetay.out","w")
print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WY", "WYERR","WYM","PHIY", "PHIYERR","PHIM"
print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le"

bpms=intersect(listy)
bpms=modelIntersect(bpms,model)
dolinregbet(filefile,dpplist,betalisty,bpms,"V","beta",zeroby,model)
filefile.close()

print "Driven beta finished"\

#
# driven coupling
#
print "Driven coupling"

filefile=open(output+"/chromcoupling.out","w")
print >>filefile,"NAME S CHROMCOUPLE  CHROMe  CHROMMDL"
print >>filefile,"%s   %le  %le       %le     %le"

bpms=intersect(couplel)
bpms=modelIntersect(bpms,model)

dolinregCoupling(couplelist,bpms,dpplist,filefile,model)
filefile.close()


print "Driven coupling finished"
filefile.close()


if freeswitch==1:
  #
  # free beta
  #
  print "Free beta"
  #H
  filefile=open(output+"/chrombetax_free.out","w")
  print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WX","WXERR","WMO","PHIX", "PHIXERR","PHIM"
  print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le"

  bpms=intersect(listx)
  bpms=modelIntersect(bpms,model)
  dolinregbet(filefile,dpplist,betalistx,bpms,"H","beta",zerobxf,model)
  filefile.close()

  #V
  filefile=open(output+"/chrombetay_free.out","w")
  print >>filefile, "* NAME", "S",  "dbb", "dbberr", "dalfa", "daerr", "WY", "WYERR","WYM","PHIY", "PHIYERR","PHIM"
  print >>filefile, "$ %s  %le  %le  %le  %le  %le %le %le %le  %le"

  bpms=intersect(listy)
  bpms=modelIntersect(bpms,model)
  dolinregbet(filefile,dpplist,betalisty,bpms,"V","beta",zerobyf,model)
  filefile.close()

  print "Free beta finished"

  #
  # free coupling
  #
  print "Free coupling"

  filefile=open(output+"/chromcoupling_free.out","w")
  print >>filefile,"NAME S CHROMCOUPLE  CHROMe  CHROMMDL"
  print >>filefile,"%s   %le  %le       %le     %le"

  bpms=intersect(couplel)
  bpms=modelIntersect(bpms,model)

  dolinregCoupling(couplelist,bpms,dpplist,filefile,model)
  filefile.close()


  print "Free coupling finished"
  filefile.close()
