################################################################
#                                                              #
#  @ Glenn Vanbavinckhove  (gvanbavi@cern.ch)=> Date: 11/09/09 #
#                                                              #
################################################################
#
#  !=> SegementBySegment_0.0.py : Construction of the program (11/09/09)
#  !=> SegementBySegment_0.1.py : - Adding the dispersion to the program (23/09/09)
#                                 - Make distuinsh between BPMs and instuments (24/09/09)
#                                 - Adding output tables (sbs...out) (21/10/09)
#                                 - Adding coupling propogation (21/10/09)
#                                 - changing output when segment is start and end (24/11/09)
#
#

###### imports
from optparse import OptionParser
from metaclass import twiss
import os,sys
from math import sqrt




###### optionparser
parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                metavar="ACCEL", default="LHCB1",dest="accel")
#parser.add_option("-d", "--dictionary",
#                help="File with the BPM dictionary",
#                metavar="DICT", default="0", dest="dict")
parser.add_option("-f", "--path", # assumes that output is same as input
                help="Path to measurement files",
                metavar="PATH", default="./", dest="path")
parser.add_option("-s", "--start",
                help="file that defines the segment",
                metavar="SEGF", default="./", dest="segf")
parser.add_option("-t", "--twiss",
                help="basic twiss file",
                metavar="TWISS", default="./", dest="twiss")
parser.add_option("-r", "--response",
                help="switch calculate corrections 0 no, 1 yes",
                metavar="COR", default="0", dest="cor")
parser.add_option("-g", "--graphic",
                help="choice between graphic or table output 0=graphic and 1=table output",
                metavar="GRA", default="1", dest="gra")
parser.add_option("-p", "--save", # assumes that output is same as input
                help="Path to measurement files",
                metavar="SAVE", default="./", dest="SAVE")
parser.add_option("-m", "--mad", # assumes that output is same as input
                help="mad link",
		  metavar="mad", default="", dest="mad")
parser.add_option("-b", "--bbsrouce", # assumes that output is same as input
                help="beta beat source",
                metavar="bb", default="./", dest="bb")
parser.add_option("-x", "--take", # assumes that output is same as input
                help="take or create madx 0/1",
                metavar="mad", default="0", dest="madpass")

(options, args) = parser.parse_args()




####### some defs
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

def getnextelement(names,mea,meay,model,end,value):

	print "Element "+element+" not found in measurement ... will look for other one!!"
	locindx=model.indx[element.upper()]
	endlocindx=model.indx[end.upper()]
	for el in names:

		llocindx=model.indx[el.upper()]
		ell="null"

		if value==0:
			if (llocindx>locindx and llocindx<endlocindx):
				try:
					if(mea.BETX[mea.indx[element.upper()]] >0 and meay.BETY[mea.indx[element.upper()]] >0):
						ell=el
						break
				except:
					print "element not in measurement"
			else:
				print "looking element ",el," is not ok"
		if value==1:
			if (llocindx>locindx and llocindx<endlocindx):
				try:
					if(mea.BETX[mea.indx[element.upper()]] >0 and meay.BETY[mea.indx[element.upper()]] >0):
						ell=el
						break
				except:
					print "element not in measurement"
			else:
				print "looking element ",el," is not ok"

		
            
        if ell=="null":

            print "No suitable candidate found ... sys exit (try to increase the boundaries)"
            sys.exit()

	return ell
	

def meascontains(mea,meay,element,end,model):

    # check first if end element is there:
    # 1 is end, 0 is beginning
    try:
	    endloc=mea.S[mea.indx[end.upper()]]
	    
    except:
	    print "end element is not in measurement !!"
	    ell=getnextelement(names,mea,meay,model,end,1)

    try:
        names=mea.NAME
        loc=mea.S[mea.indx[element.upper()]]
        endloc=mea.S[mea.indx[end.upper()]]
    
        check=mea.indx[element.upper()]
	if(mea.BETX[mea.indx[element.upper()]] >0 and meay.BETY[mea.indx[element.upper()]] >0):
	
		ell=element
	else:
		print "beta at bpm ",element," is containing a negative beta ... will check for other element"
		ell=getnextelement(names,mea,meay,model,end)
    except:

        ell=getnextelement(names,mea,meay,model,end)

    return ell

def getelement(filedatax,filedatay,eelement,model):

    loc=model.S[model.indx[eelement]]
    
    names=filedatax.NAME

    pos=-100
    epos=3000000
    name="null"
    ename="null"
    
    print "Looking for element just in front of : ",eelement, " at location : ",str(loc)

    for el in names:

	    local=model.S[model.indx[el]]

	    try:

		    if(filedatax.BETX[filedatax.indx[el.upper()]] >0 and filedatay.BETY[filedatay.indx[el.upper()]] >0):
			    if local<loc and local>pos:

				    pos=local
				    name=el
			    
			    if local>loc and local<epos:
				    epos=local
				    ename=el
	    except:
		    print "bpm ",el," not found in both planes"

    if pos==-100 or epos==3000000:
	    print "no element found !!! => system exit"
	    sys.exit()
			    
    print "element infront found at ",str(pos)," with name ",name," at a distance of ", str(loc-pos)," from ",eelement
    print "element after found at ",str(epos)," with name ",ename," at a distance of ", str(epos-loc)," from ",eelement
    
    return name,ename

def getValues(filedatax,filedatay,element,dp,filedatadx,filedatady):

    betax=filedatax.BETX[filedatax.indx[element]]
    errbetax=filedatax.ERRBETX[filedatax.indx[element]]
    alfx=filedatax.ALFX[filedatax.indx[element]]
    erralfx=filedatax.ERRALFX[filedatax.indx[element]]
    
    betay=filedatay.BETY[filedatay.indx[element]]
    errbetay=filedatay.ERRBETY[filedatay.indx[element]]
    alfy=filedatay.ALFY[filedatay.indx[element]]
    erralfy=filedatay.ERRALFY[filedatay.indx[element]]

    hor=[betax,errbetax,alfx,erralfx]
    ver=[betay,errbetay,alfy,erralfy]


    
    if (dp==0):

	    dp=[0,0,0,0,0,0]

    elif (dp==1):

	    try:
		    dx=filedatadx.DX[filedatadx.indx[element]]
		    dpx=filedatadx.DPX[filedatadx.indx[element]]
		    dxeerx=filedatadx.STDDX[filedatadx.indx[element]]
	    except:
		    dx=0
		    dpx=0
		    dxeerx=0
	    try:
		    dy=filedatady.DY[filedatady.indx[element]]
		    dyeery=filedatady.STDDY[filedatady.indx[element]]	    
		    dpy=filedatady.DPY[filedatady.indx[element]]
	    except:
		    dy=0
		    dyeery=0
		    dpy=0
	    
	    #dy=0
	    #dpy=0

	    dp=[dx,dpx,dy,dpy,dxeerx,dyeery]

    return[hor,ver,dp]

def getTwiss(filee,element):

    filedatax=twiss(filee)

    betax=filedatax.BETX[filedatax.indx[element]]
    alfx=filedatax.ALFX[filedatax.indx[element]]
    sx=filedatax.S[filedatax.indx[element]]

    betay=filedatax.BETY[filedatax.indx[element]]
    alfy=filedatax.ALFY[filedatax.indx[element]]
    sy=filedatax.S[filedatax.indx[element]]

    dx=filedatax.DX[filedatax.indx[element]]
    dpx=filedatax.DPX[filedatax.indx[element]]
    dy=filedatax.DY[filedatax.indx[element]]
    dpy=filedatax.DPY[filedatax.indx[element]]
 
    hor=[betax,alfx,sx]
    ver=[betay,alfy,sy]
    dp=[dx,dpx,dy,dpy]

    return[hor,ver,dp]

def run4mad(path,hor,ver,hore,vere,dp,dpe):

    if options.accel=="LHCB2":

        dire=-1
        start="MKI.A5R8.B2"
        beam="B2"
        
    elif options.accel=="LHCB1":

        dire=1
        start="MKI.A5L2.B1"
        beam="B1"

    cpath=options.bb
     
    filename=path+'/var4mad.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%BETX/\''+str(hor[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETY/\''+str(ver[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETX/\''+str(hor[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETY/\''+str(ver[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFX/\''+str(hor[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFY/\''+str(ver[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFX/\''+str(hor[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFY/\''+str(ver[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DX/\''+str(dp[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDX/\''+str(dp[4])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DY/\''+str(dp[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDY/\''+str(dp[5])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPX/\''+str(dp[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPY/\''+str(dp[3])+'\'/g\' \\\n')
    
    file4nad.write('    -e \'s/%ENDBX/\''+str(hore[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDBY/\''+str(vere[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRENDBX/\''+str(hore[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRENDBY/\''+str(vere[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFENDX/\''+str(-hore[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ALFENDY/\''+str(-vere[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFENDX/\''+str(hore[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRALFENDY/\''+str(vere[3])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DENDX/\''+str(dpe[0])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDENDX/\''+str(dpe[4])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DENDY/\''+str(dpe[2])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRDENDY/\''+str(dpe[5])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPENDX/\''+str(-dpe[1])+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DPENDY/\''+str(-dpe[3])+'\'/g\' \\\n')
    
    file4nad.write('    -e \'s/%STARTFROM/\''+str(element)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDAT/\''+str(eelement)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(file4seg.LABEL)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DIRE/\''+str(dire)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%START/\''+str(start)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BEAM/\''+str(beam)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PATH/\'\"'+str(path.replace('/','\/'))+'\"\'/g\' \\\n')
    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/job.InterpolateBetas.mask > '+path+'/t_'+str(file4seg.LABEL)+'.madx \n')

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))

    runmad(path)
   
    
def runmad(path):
	
	os.system(options.mad+'madx < '+path+'t_'+str(file4seg.LABEL)+'.madx')
	
   
  
def run4plot(path,spos,epos,beta4plot,cpath,meapath):

    filename=path+'/var4plot.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%PATH/\'\"'+str(path.replace("/","\/"))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%EndPoint/\''+str(epos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%StartPoint/\''+str(spos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(file4seg.LABEL)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETA/\''+str(beta4plot)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%MEA/\'\"'+str(meapath.replace("/","\/"))+'\"\'/g\' \\\n')    
    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.mask > '+path+'/gplot_'+file4seg.LABEL)

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))
    os.system("chmod 777 "+str(path+'/gplot_'+file4seg.LABEL))
    os.system(str(path+'/gplot_'+file4seg.LABEL))

 
   
def GetPhaseEM(exp, mod):
    phasem=[]
    bpm1=[]
    bpm2=[]
    s1=[]
    s2=[]
    phaseexp=[]
    for elm in mod.NAME:
        status=1
        try:
            expind=exp.indx[elm]
        except:
            print elm, "not found in exp"
            status=0

        if status==1:
            el2=exp.NAME2[exp.indx[elm]]
            elmind=mod.indx[elm]
            try:
                elm2ind=mod.indx[el2]
            except:
                print el2, "not found in model"
                status=0

            if status==1:
                if "PHXMDL" in  exp.__dict__.keys():
                    modphaseadv=mod.MUX[elm2ind]-mod.MUX[elmind]
		    phaseexp.append(exp.PHASEX[expind])
                    
                elif "PHYMDL" in  exp.__dict__.keys():
                    modphaseadv=mod.MUY[elm2ind]-mod.MUY[elmind]
		    phaseexp.append(exp.PHASEY[expind])
                    
                    
                bpm1.append(elm)
                bpm2.append(mod.NAME[elm2ind])
                s1.append(mod.S[elmind])
                s2.append(mod.S[elm2ind])
                phasem.append(modphaseadv)
                
                
    return bpm1, bpm2, s1, s2, phaseexp, phasem

def writePhase(filename,bpm1, bpm2, s1, s2, phaseexp, phasem ):

    
    f=open(filename, "w")
    for i in range(len(bpm1)):
        print >>f, bpm1[i], bpm2[i], s1[i], s2[i], phaseexp[i], phasem[i]

    f.close()


def getcoupling(couplingtwiss,model,beginelement,endelement):

	modell=twiss(model)
	modell.Cmatrix()

	F1001=[]
	F1010=[]

	#begin
	f1001ab=couplingtwiss.F1001W[couplingtwiss.indx[beginelement]]
	f1001rb=couplingtwiss.F1001R[couplingtwiss.indx[beginelement]]
	f1001ib=couplingtwiss.F1001I[couplingtwiss.indx[beginelement]]

	f1001abm=abs(modell.f1001[modell.indx[beginelement]])
	f1001rbm=modell.f1001[modell.indx[beginelement]].real
	f1001ibm=modell.f1001[modell.indx[beginelement]].imag

	f1010ab=couplingtwiss.F1010W[couplingtwiss.indx[beginelement]]
	f1010rb=couplingtwiss.F1010R[couplingtwiss.indx[beginelement]]
	f1010ib=couplingtwiss.F1010I[couplingtwiss.indx[beginelement]]

	f1010abm=abs(modell.f1010[modell.indx[beginelement]])
	f1010rbm=modell.f1010[modell.indx[beginelement]].real
	f1010ibm=modell.f1010[modell.indx[beginelement]].imag

	#end
	f1001ae=couplingtwiss.F1001W[couplingtwiss.indx[endelement]]
	f1001re=couplingtwiss.F1001R[couplingtwiss.indx[endelement]]
	f1001ie=couplingtwiss.F1001I[couplingtwiss.indx[endelement]]

	f1001aem=abs(modell.f1001[modell.indx[endelement]])
	f1001rem=modell.f1001[modell.indx[endelement]].real
	f1001iem=modell.f1001[modell.indx[endelement]].imag

	f1010ae=couplingtwiss.F1010W[couplingtwiss.indx[endelement]]
	f1010re=couplingtwiss.F1010R[couplingtwiss.indx[endelement]]
	f1010ie=couplingtwiss.F1010I[couplingtwiss.indx[endelement]]

	f1010aem=abs(modell.f1010[modell.indx[endelement]])
	f1010rem=modell.f1010[modell.indx[endelement]].real
	f1010iem=modell.f1010[modell.indx[endelement]].imag	
	
	F1001=[f1001ab,f1001rb,f1001ib,f1001abm,f1001rbm,f1001ibm,f1001ae,f1001re,f1001ie,f1001aem,f1001rem,f1001iem]

	F1010=[f1010ab,f1010rb,f1010ib,f1010abm,f1010rbm,f1010ibm,f1010ae,f1010re,f1010ie,f1010aem,f1010rem,f1010iem]

	return [F1001,F1010]
	

def reversetable(path):
	file=open(path+"/twiss_"+file4seg.LABEL+"_back_rev.dat",'w')
	base=twiss(path+"/twiss_"+file4seg.LABEL+"_back.dat")

	file.write("* NAME                                S               BETX               ALFX               BETY               ALFY    DX       DPX     DY     DPY\n")
	file.write("$ %s                                %le                %le                %le                %le                %le      %le                %le                %le                %le \n")
	
	bpms=base.NAME
	s=base.S
	endpos=base.S[base.indx[bpms[len(bpms)-1]]]
	betx=base.BETX
	bety=base.BETY
	alfx=base.ALFX
	alfy=base.ALFY
	dx=base.DX
	dpx=base.DPX
	dy=base.DY
	dpy=base.DPY

	for i in range(len(bpms)):
		index=len(bpms)-1-i
		bpm=bpms[index]
		bex=betx[index]
		bey=bety[index]
		alx=alfx[index]
		aly=alfy[index]
		dex=dx[index]
		depx=dpx[index]
		dey=dy[index]
		depy=dpy[index]
		ss=s[index]

		print bpm, dex

		file.write(str(bpm)+' '+str(endpos-ss)+' '+str(bex)+' '+str(alx)+' '+str(bey)+' '+str(aly)+'  '+str(dex)+' '+str(depx)+' '+str(dey)+' '+str(depy)+'\n')
		
	file.close()


def rotateparts(bpms,img1001,real1001,img1010,real1010,phasex,phasey):

	img1001l=[]
	real1001l=[]
	f1001=[]
	p1001=[]
	img1010l=[]
	real1010l=[]
	f1010=[]
	p1010=[]

	imgl=[]
	reall=[]
	
	img1001l.append(img1001)
	real1001l.append(real1001)
	img1010l.append(img1010)
	real1010l.append(real1010)
	
	
	for bpm in range(len(bpms)):

		if bpm<len(bpms):

			namebpm=bpms[bpm+1].upper()

			phx=phasex.PHASEX[phasex.indx[namebpm]]
			phy=phasey.PHASEY[phasey.indx[namebpm]]
			delta1001=(phx-phy)
			delta1010=(phx+phy)
		
			img1001t=img1001*complex(cos(delta1001),-sin(delta1001))
			img1001l.append(img1001t)
			real1001t=real1001*complex(cos(delta1001),-sin(delta1001))
			real1001l.append(real1001t)
			img1010t=img1010*complex(cos(delta1010),-sin(delta1010))
			img1010l.append(img1010t)
			real1010t=real1010*complex(cos(delta1010),-sin(delta1010))
			real1010l.append(real1010t)
		
			f1001.append(sqrt(img1001t**2+real1001t**2))
			p1001.append(atan(img1001t/real1001t)%1)

			f1010.append(sqrt(img1010t**2+real1010t**2))
			p1010.append(atan(img1010t/real1010t)%1)

		

	imgl=[img1001l,img1010l]
	reall=[real1001l,real1010l]

	return [bpms,imgl,reall]


def createTables(outputname,path,columnnames,paranames,data,mainvariable,mainvalue):

	filefile=open(path+"/"+outputname,'w')

	# writing main variables
	if len(mainvariable)==len(mainvalue):
		for count in range(len(mainvariable)):
			
			if isinstance(mainvalue[count],str):
				filefile.write('@ '+mainvariable[count]+' %s '+str(mainvalue[count])+'\n')
			
			else:
				filefile.write('@ '+mainvariable[count]+' %le '+str(mainvalue[count])+'\n')
	else:
		print "cannot write main variables ..."

	# writing columnnames and paranames

	if len(columnnames)==len(paranames):


		filefile.write('  '.join(columnnames)+'\n')
		filefile.write('  '.join(paranames)+'\n')

		#writing data
		#
		# data[x][y]
		#
		#   x = data set
		#   y = data in data set
		#   => has to bed added in same order
		#
		
		for y in range(len(data)):
		
			filefile.write(data[y]+'\n')


	
	
	else:

		print "cannot create table names for columns are not equal => system exit "
		sys.exit()

	filefile.close()

####### main part
file4seg=twiss(options.segf)
path=options.path
filedatax=twiss(path+"/getbetax.out")
betxtwiss=filedatax
filedatay=twiss(path+"/getbetay.out")
betytwiss=filedatay
filedatax=twiss(path+"/getampbetax.out")
ampbetxtwiss=filedata
filedatay=twiss(path+"/getampbetay.out")
ampbetytwiss=filedatay
filephasex=twiss(path+"/getphasey.out")
filephasey=twiss(path+"/getphasex.out")
#filecoupling=twiss(path+"getcoupling.out")
try:
	filedx=twiss(path+"/getDx.out")
	filedy=twiss(path+"/getDy.out")
	disp=1
except:
	print "no dispersion files... will continue without taking into account dispersion"
	filedx=[]
	filedy=[]
	disp=0





# checking if end element is BPM or instrument (IP,collimators)
list2run=[]
elementswitch=-1
print options.segf.split('.tfs')[0]
if "COLLI" in options.segf.split('.tfs')[0]:
#	sys.exit()
	list2run=file4seg.NAME
	if len(list2run)>0:
		if "T" in list2run[0]:
			print "end element is collimator : ",list2run[0]
			elementswitch=1
#		elif "IP" in names[0]:
#			print "end element is interation region : ",list2run[0]
#			elementswitch=2
	else:
		print "list doesnt contain any input, please check => ",options.segf
		sys.exit()
elif "IP" in options.segf.split('.tfs')[0]:
	elementswitch=0
	print "end element is interation region : ",file4seg.EBPM
	list2run.append(file4seg.EBPM)
	
elif "BPM" in file4seg.EBPM:
	print "end element is BPM : ",file4seg.EBPM
	list2run.append(file4seg.EBPM)
	elementswitch=0
	print " Calculating the twiss functions from "+file4seg.SBPM+" until "+file4seg.EBPM
else:
	print "element is neither BPM, collimator or IP :",file4seg.EBPM
	print "Please check => system exit"
	sys.exit()
mainvariablebetx=['NFILES']
mainvaluebetx=[filedatax.COUNT[0]]
columnnamesbetx=['* NAME','S','BETX','ERRBETX','BETXMDL','BETXP','ERRBETXP','BETXB','ERRBETXB']
variablenamebetx=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
databetx=[]
# for ampbetx
mainvariableampbetx=['NFILES']
mainvalueampbetx=[filedatax.COUNT[0]]
columnnamesampbetx=['* NAME','S','BETX','ERRBETX','BETXMDL','BETXP','ERRBETXP','BETXB','ERRBETXB']
variablenameampbetx=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
dataampbetx=[]
# for alfx
columnnamesalfx=['* NAME','S','ALFX','ERRALFX','ALFXMDL','ALFXP','ERRALFXP','ALFXB','ERRALFXB']
variablenamealfx=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablealfx=['NFILES']
mainvaluealfx=[filedatax.COUNT[0]]
dataalfx=[]
# for bety
betytwiss=twiss(path+"/getbetay.out")
columnnamesbety=['* NAME','S','BETY','ERRBETY','BETYMDL','BETYP','ERRBETYP','BETYB','ERRBETYB']
variablenamebety=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablebety=['NFILES']
mainvaluebety=[filedatay.COUNT[0]]
databety=[]
# for ampbety
columnnamesampbety=['* NAME','S','BETY','ERRBETY','BETYMDL','BETYP','ERRBETYP','BETYB','ERRBETYB']
variablenameampbety=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariableampbety=['NFILES']
mainvalueampbety=[filedatay.COUNT[0]]
dataampbety=[]
# for alfy
columnnamesalfy=['* NAME','S','ALFY','ERRALFY','ALFYMDL','ALFYP','ERRALFYP','ALFYB','ERRALFYB']
variablenamealfy=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablealfy=['NFILES']
mainvaluealfy=[filedatay.COUNT[0]]
dataalfy=[]
#for dx and dpx
dx4twiss=path+"/getDx.out"
columnnames1dx=['* NAME','S','DX','ERRDX','DXMDL','DXP','ERRDXP','DXB','ERRDXB']
variablename1dx=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
columnnames2dx=['* NAME','S','DPX','ERRDPX','DPXMDL','DPXP','ERRDPXP','DPXB','ERRDPXB']
variablename2dx=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariabledx=['NFILES']
mainvaluedx=[0]
datadx=[]
datadpx=[]
# for dy and dpy
dy4twiss=path+"/getDy.out"
columnnames1dy=['* NAME','S','DY','ERRDY','DYMDL','DYP','ERRDYP','DYB','ERRDYB']
variablename1dy=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
columnnames2dy=['* NAME','S','DPY','ERRDPY','DPYMDL','DPYP','ERRDPYP','DPXB','ERRDPYB']
variablename2dy=['$ %s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariabledy=['NFILES']
mainvaluedy=[0]
datady=[]
datadpy=[]
twisstwiss=twiss(options.twiss)

######## big loop

for namename in list2run:

	if elementswitch==0:
		element=file4seg.SBPM
		eelement=file4seg.EBPM
		label=file4seg.LABEL
		element=meascontains(filedatax,filedatay,element,eelement,twisstwiss)
		print "ddd "+element
		#sys.exit()
	#elif elementswitch==2:
	#	element=file4seg.SBPM
		#eelement=file4seg.EBPM
	#	label=file4seg.LABEL
	#	element=meascontains(filedatax,filedatay,element,eelement,twiss(options.twiss))
	else:
		label=file4seg.LABEL
		element,eelement=getelement(filedatax,filedatay,namename,twisstwiss)


	[hor,ver,dp]=getValues(filedatax,filedatay,element,disp,filedx,filedy)
	[hore,vere,dpe]=getValues(filedatax,filedatay,eelement,disp,filedx,filedy)

         #[F1001,F1010]=getcoupling(filecoupling,options.twiss,element,eelement)
	savepath=options.SAVE
	print options.madpass
	#sys.exit()
	if str(options.madpass)=="0":
		run4mad(savepath,hor,ver,hore,vere,dp,dpe)
	else:
		runmad(savepath)
		print "skipping mad"
	

	reversetable(savepath)

	##################
	#=> switch only for element - to - element
	if elementswitch==0:
	##################
		[hor,ver,dp]=getTwiss(savepath+"/StartPoint.twiss",element)
		startpos=hor[2]
		[hor,ver,dp]=getTwiss(savepath+"/twiss_"+str(file4seg.LABEL)+".dat",element)
		betx=hor[0]
		bety=ver[0]
		ax=hor[1]
		ay=ver[1]
		dx=dp[0]
		dpx=dp[1]
		dy=dp[2]
		dpy=dp[3]
		[hor,ver,dp]=getTwiss(savepath+"/EndPoint.twiss",eelement)
		endpos=hor[2]

		[hor,ver,dp]=getTwiss(savepath+"/twiss.b+.dat",element)
		bplusx=hor[0]
		bplusy=ver[0]
		[hor,ver,dp]=getTwiss(savepath+"/twiss.b-.dat",element)
		bminx=hor[0]
		bminy=ver[0]

		[hor,ver,dp]=getTwiss(savepath+"/twiss.a+.dat",element)
		aplusx=hor[1]
		aplusy=ver[1]
		[hor,ver,dp]=getTwiss(savepath+"/twiss.a-.dat",element)
		aminx=hor[1]
		aminy=ver[1]

		berrx=bplusx-bminx
	        aerrx=aplusx-aminx
		berry=bplusy-bminy
		aerry=aplusy-aminy

		phx=twiss(path+'/getphasex.out')
		phy=twiss(path+'/getphasey.out')

		m=twiss(savepath+"/twiss_"+file4seg.LABEL+".dat")

		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phx, m)
		writePhase(savepath+"/phasexEM.out",bpm1, bpm2, s1, s2, phaseexp, phasem)
		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phy, m)
		writePhase(savepath+"/phaseyEM.out",bpm1, bpm2, s1, s2, phaseexp, phasem)

		m=twiss(savepath+"/twiss_"+file4seg.LABEL+"_play.dat")

		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phx, m)
		writePhase(savepath+"/phasexEM_play.out",bpm1, bpm2, s1, s2, phaseexp, phasem)
		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phy, m)
		writePhase(savepath+"/phaseyEM_play.out",bpm1, bpm2, s1, s2, phaseexp, phasem)


		m=options.twiss
		[hor,ver,dp]=getTwiss(m,element)
		beta4plot=hor[2]

	        ######################
		if options.gra=='1':
			run4plot(savepath,startpos,endpos,beta4plot,options.bb,path)


		################
	        fileresul=open(savepath+'/resul_'+file4seg.LABEL+'.tfs','w')
	        fileresul.write('* NAME   S     BETX    ERRBETX    ALFX   ERRALFX   BETY   ERRBETY  ALFY   ERRALFX  DX  DPX   DY    DPY\n')
	        fileresul.write('* %s     %le    %le   %le         %le    %le       %le    %le      %le    %le     %le   %le    %le    %le\n')
	        fileresul.write(str(file4seg.EBPM)+' '+str(endpos)+' '+str(betx)+' '+str(berrx)+' '+str(ax)+' '+str(aerrx)+' '+str(bety)+' '+str(berry)+' '+str(ay)+' '+str(aerry)+' '+str(dx)+' '+str(dpx)+' '+str(dy)+' '+str(dpy))
	        fileresul.close()


	        #### creating tables
		#
		# => opening all the results
		#


		# results from propagation
		normal_pro=twiss(savepath+'/twiss_'+label+'.dat')
		back_pro=twiss(savepath+'/twiss_'+label+'_back_rev.dat')
		bpmsbetx=intersect([filedatax,normal_pro,back_pro])
		ampbpmsbetx=intersect([ampbetxtwiss,normal_pro,back_pro])
		bpmsbety=intersect([filedatay,normal_pro,back_pro])
		ampbpmsbety=intersect([ampbetytwiss,normal_pro,back_pro])
		if os.path.exists(dx4twiss):
			try:
				Dx=twiss(path+"/getDx.out")
				bpmsdx=intersect([Dx,normal_pro,back_pro])
				mainvaluedx=[Dx.COUNT[0]]
			except:
				print "no dispersion file"
				bpmsdx=[]
				mainvaluedx=[0]
		if os.path.exists(dy4twiss):
			try:
				Dy=twiss(path+"/getDy.out")
				bpmsdy=intersect([Dy,normal_pro,back_pro])
				mainvaluedy=[Dy.COUNT[0]]
			except:
				print "no dispersion file"
				bpmsdy=[]
				mainvaluedy=[0]
		
		errbetamin=twiss(savepath+'/twiss.b-.dat')
		errbetaminb=twiss(savepath+'/twiss.b-_back.dat')
	
		errbetamax=twiss(savepath+'/twiss.b+.dat')
		errbetamaxb=twiss(savepath+'/twiss.b+_back.dat')

		erralfmin=twiss(savepath+'/twiss.a-.dat')
		erralfminb=twiss(savepath+'/twiss.a-_back.dat')

		erralfmax=twiss(savepath+'/twiss.a+.dat')
		erralfmaxb=twiss(savepath+'/twiss.a+_back.dat')

		errdmin=twiss(savepath+'/twiss.d-.dat')
		errdminb=twiss(savepath+'/twiss.d-_back.dat')

		errdmax=twiss(savepath+'/twiss.d+.dat')
		errdmaxb=twiss(savepath+'/twiss.d+_back.dat')


		# results from getllm
                ##
		#
		#betx
		print "filling table for betx"
		for bpm in bpmsbetx:
			name=bpm[1]
			s=betxtwiss.S[betxtwiss.indx[name]]
			bet=betxtwiss.BETX[betxtwiss.indx[name]]
			errbet=sqrt(betxtwiss.ERRBETX[betxtwiss.indx[name]]**2+betxtwiss.STDBETX[betxtwiss.indx[name]]**2)
			betmdl=betxtwiss.BETXMDL[betxtwiss.indx[name]]

			betp=normal_pro.BETX[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETX[errbetamax.indx[name]]-errbetamin.BETX[errbetamin.indx[name]])

			betb=back_pro.BETX[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETX[errbetamaxb.indx[name]]-errbetaminb.BETX[errbetaminb.indx[name]])
	
			databetx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		# for ampbetx
			print "filling table for betx"
		for bpm in ampbpmsbetx:
			name=bpm[1]
			s=ampbetxtwiss.S[ampbetxtwiss.indx[name]]
			bet=ampbetxtwiss.BETX[ampbetxtwiss.indx[name]]
			errbet=ampbetxtwiss.STDBETX[ampbetxtwiss.indx[name]]
			betmdl=ampbetxtwiss.BETXMDL[ampbetxtwiss.indx[name]]

			betp=normal_pro.BETX[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETX[errbetamax.indx[name]]-errbetamin.BETX[errbetamin.indx[name]])

			betb=back_pro.BETX[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETX[errbetamaxb.indx[name]]-errbetaminb.BETX[errbetaminb.indx[name]])
	
			dataampbetx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		
		# for alphax
		print "filling table for alfx"
		for bpm in bpmsbetx:
			name=bpm[1]
			s=betxtwiss.S[betxtwiss.indx[name]]
			bet=betxtwiss.ALFX[betxtwiss.indx[name]]
			errbet=sqrt(betxtwiss.ERRALFX[betxtwiss.indx[name]]**2+betxtwiss.STDALFX[betxtwiss.indx[name]]**2)
			betmdl=betxtwiss.ALFXMDL[betxtwiss.indx[name]]
			
			betp=normal_pro.ALFX[normal_pro.indx[name]]
			errbetp=abs(erralfmax.ALFX[erralfmax.indx[name]]-erralfmin.ALFX[erralfmin.indx[name]])

			betb=-back_pro.ALFX[back_pro.indx[name]]
			errbetb=abs(erralfmaxb.ALFX[erralfmaxb.indx[name]]-erralfminb.ALFX[erralfminb.indx[name]])
	
			dataalfx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
					
                #for bety
		print "filling table for bety"
		for bpm in bpmsbety:
			name=bpm[1]
			s=betytwiss.S[betytwiss.indx[name]]
			bet=betytwiss.BETY[betytwiss.indx[name]]
			errbet=sqrt(betytwiss.ERRBETY[betytwiss.indx[name]]**2+betytwiss.STDBETY[betytwiss.indx[name]]**2)
			betmdl=betytwiss.BETYMDL[betytwiss.indx[name]]

			betp=normal_pro.BETY[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETY[errbetamax.indx[name]]-errbetamin.BETY[errbetamin.indx[name]])

			betb=back_pro.BETY[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETY[errbetamaxb.indx[name]]-errbetaminb.BETY[errbetaminb.indx[name]])
	
			databety.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
			
		# for ampbety
		for bpm in ampbpmsbety:
			name=bpm[1]
			s=ampbetytwiss.S[ampbetytwiss.indx[name]]
			bet=ampbetytwiss.BETY[ampbetytwiss.indx[name]]
			errbet=ampbetytwiss.STDBETY[ampbetytwiss.indx[name]]
			betmdl=betytwiss.BETYMDL[ampbetytwiss.indx[name]]

			betp=normal_pro.BETY[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETY[errbetamax.indx[name]]-errbetamin.BETY[errbetamin.indx[name]])

			betb=back_pro.BETY[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETY[errbetamaxb.indx[name]]-errbetaminb.BETY[errbetaminb.indx[name]])
	
			databety.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		
		# for alfy
		print "filling table for alfy"
		for bpm in bpmsbety:
			name=bpm[1]
			s=betytwiss.S[betytwiss.indx[name]]
			bet=betytwiss.ALFY[betytwiss.indx[name]]
			errbet=sqrt(betytwiss.ERRALFY[betytwiss.indx[name]]**2+betytwiss.STDALFY[betytwiss.indx[name]]**2)
			betmdl=betytwiss.ALFYMDL[betytwiss.indx[name]]

			betp=normal_pro.ALFY[normal_pro.indx[name]]
			errbetp=abs(erralfmax.ALFY[erralfmax.indx[name]]-erralfmin.ALFY[erralfmin.indx[name]])
			
			betb=-back_pro.ALFY[back_pro.indx[name]]
			errbetb=abs(erralfmaxb.ALFY[erralfmaxb.indx[name]]-erralfminb.ALFY[erralfminb.indx[name]])
			dataalfy.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		
	         #for Dx and Dpx
		print "filling table for dx and dpx"
		for bpm in bpmsdx:
			name=bpm[1]
			s=Dx.S[Dx.indx[name]]
			dx=Dx.DX[Dx.indx[name]]
		        errdx=Dx.STDDX[Dx.indx[name]]
		        dxmdl=Dx.DXMDL[Dx.indx[name]]
			 
		        dpx=Dx.DPX[Dx.indx[name]]
		        dpxmdl=Dx.DPXMDL[Dx.indx[name]]
		        dpxerr="0"
			 
		        dxp=normal_pro.DPX[normal_pro.indx[name]]
		        errdxp=abs(errdmax.DX[errdmax.indx[name]]-errdmin.DX[errdmin.indx[name]])
		
		        dpxp=normal_pro.DPX[normal_pro.indx[name]]
		        dpxperr="0"
			 
		        dxb=back_pro.DX[back_pro.indx[name]]
		        errdxb=abs(errdmaxb.DX[errdmaxb.indx[name]]-errdminb.DX[errdminb.indx[name]])
		
		        dpxb=-back_pro.DPX[back_pro.indx[name]]
		        dpxberr="0"
	
			datadx.append(name+" "+str(s)+" "+str(dx)+" "+str(errdx)+" "+str(dxmdl)+" "+str(dxp)+" "+str(errdxp)+" "+str(dxb)+" "+str(errdxb))
			datadpx.append(name+" "+str(s)+" "+str(dpx)+" "+str(dpxerr)+" "+str(dpxmdl)+" "+str(dpxp)+" "+str(dpxperr)+" "+str(dpxb)+" "+str(dpxberr))

		 # for Dy and Dpy
		print "filing table for dy and dpy"
	        for bpm in bpmsdy:
			name=bpm[1]
		        s=Dy.S[Dy.indx[name]]
		        dx=Dy.DY[Dy.indx[name]]
		        errdx=Dy.STDDY[Dy.indx[name]]
		        dxmdl=Dy.DYMDL[Dy.indx[name]]
			 
		        dpx=Dy.DPY[Dy.indx[name]]
		        dpxmdl=Dy.DPYMDL[Dy.indx[name]]
		        dpxerr="0"
		        
		        dxp=normal_pro.DPY[normal_pro.indx[name]]
		        errdxp=abs(errdmax.DY[errdmax.indx[name]]-errdmin.DY[errdmin.indx[name]])

		        dpxp=normal_pro.DPY[normal_pro.indx[name]]
		        dpxperr="0"
			 
		        dxb=back_pro.DY[back_pro.indx[name]]
		        errdxb=abs(errdmaxb.DY[errdmaxb.indx[name]]-errdminb.DY[errdminb.indx[name]])

		        dpxb=back_pro.DPY[-back_pro.indx[name]]
		        dpxberr="0"

			datady.append(name+" "+str(s)+" "+str(dx)+" "+str(errdx)+" "+str(dxmdl)+" "+str(dxp)+" "+str(errdxp)+" "+str(dxb)+" "+str(errdxb))
			datadpy.append(name+" "+str(s)+" "+str(dpx)+" "+str(dpxerr)+" "+str(dpxmdl)+" "+str(dpxp)+" "+str(dpxperr)+" "+str(dpxb)+" "+str(dpxberr))

                 #########

		


                 #=> switch only for element - to - instrument
	elif elementswitch!=0:

		name=namename

		errbetamin=twiss(savepath+'/twiss.b-.dat')
		errbetaminb=twiss(savepath+'/twiss.b-_back.dat')
	
		errbetamax=twiss(savepath+'/twiss.b+.dat')
		errbetamaxb=twiss(savepath+'/twiss.b+_back.dat')

		erralfmin=twiss(savepath+'/twiss.a-.dat')
		erralfminb=twiss(savepath+'/twiss.a-_back.dat')

		erralfmax=twiss(savepath+'/twiss.a+.dat')
		erralfmaxb=twiss(savepath+'/twiss.a+_back.dat')

		errdmin=twiss(savepath+'/twiss.d-.dat')
		errdminb=twiss(savepath+'/twiss.d-_back.dat')

		errdmax=twiss(savepath+'/twiss.d+.dat')
		errdmaxb=twiss(savepath+'/twiss.d+_back.dat')

	

		fileresul=open(savepath+'/resul_'+namename+'.tfs','w')
		[hor,ver,dp]=getTwiss(savepath+"/twiss_"+str(file4seg.LABEL)+".dat",namename)
		[horb,verb,dpb]=getTwiss(savepath+"/twiss_"+file4seg.LABEL+"_back_rev.dat",namename)
		print "print the values"


		fileresul.write("@ INSTRUMENT  %s "+namename+"\n")
		fileresul.write("@ S %le  "+str(hor[2])+"\n")

	
		fileresul.write("* METHOD  BETX  ALFX   BETY   ALFY   DX   DPX   DY   DPY\n")
		fileresul.write("$  %s  %le  %le  %le  %le  %le  %le  %le  %le\n")

		fileresul.write("normal "+str(hor[0])+" "+str(hor[1])+" "+str(ver[0])+" "+str(ver[1])+" "+str(dp[0])+" "+str(dp[1])+" "+str(dp[2])+" "+str(dp[3])+"\n")
		fileresul.write("back "+str(horb[0])+" "+str(horb[1])+" "+str(verb[0])+" "+str(verb[1])+" "+str(dpb[0])+" "+str(dpb[1])+" "+str(dpb[2])+" "+str(dpb[3])+"\n")
		fileresul.close()

			        #### creating tables
		#
		# => opening all the results
		#

		# opening data
		normal_pro=twiss(savepath+'/twiss_'+label+'.dat')
		back_pro=twiss(savepath+'/twiss_'+label+'_back_rev.dat')

		# for betx
		s=twisstwiss.S[twisstwiss.indx[name]]
		bet=0
		errbet=0
		betmdl=twisstwiss.BETX[twisstwiss.indx[name]]
		betp=normal_pro.BETX[normal_pro.indx[name]]
		errbetp=abs(errbetamax.BETX[errbetamax.indx[name]]-errbetamin.BETX[errbetamin.indx[name]])
		betb=back_pro.BETX[back_pro.indx[name]]
		errbetb=abs(errbetamaxb.BETX[errbetamaxb.indx[name]]-errbetaminb.BETX[errbetaminb.indx[name]])
		databetx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))

	
		# for bety
		s=twisstwiss.S[twisstwiss.indx[name]]
		bet=0
		errbet=0
		betmdl=twisstwiss.BETY[twisstwiss.indx[name]]
		betp=normal_pro.BETY[normal_pro.indx[name]]
		errbetp=abs(errbetamax.BETY[errbetamax.indx[name]]-errbetamin.BETY[errbetamin.indx[name]])
		betb=back_pro.BETY[back_pro.indx[name]]
		errbetb=abs(errbetamaxb.BETY[errbetamaxb.indx[name]]-errbetaminb.BETY[errbetaminb.indx[name]])
		databety.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))

		# for alfx
		s=twisstwiss.S[twisstwiss.indx[name]]
		bet=0
		errbet=0
		betmdl=twisstwiss.ALFX[twisstwiss.indx[name]]
		betp=normal_pro.ALFX[normal_pro.indx[name]]
		errbetp=abs(erralfmax.ALFX[erralfmax.indx[name]]-erralfmin.ALFX[erralfmin.indx[name]])
		betb=-back_pro.ALFX[back_pro.indx[name]]
		errbetb=abs(erralfmaxb.ALFX[erralfmaxb.indx[name]]-erralfminb.ALFX[erralfminb.indx[name]])
		dataalfx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		# for alfy
		betmdl=twisstwiss.ALFY[twisstwiss.indx[name]]
		betp=normal_pro.ALFY[normal_pro.indx[name]]
		errbetp=abs(erralfmax.ALFY[erralfmax.indx[name]]-erralfmin.ALFY[erralfmin.indx[name]])
	      	betb=-back_pro.ALFY[back_pro.indx[name]]
		errbetb=abs(erralfmaxb.ALFY[erralfmaxb.indx[name]]-erralfminb.ALFY[erralfminb.indx[name]])
		dataalfy.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		# for dx and dpx
		s=twisstwiss.S[twisstwiss.indx[name]]
		dx=0
		dpx=0
		dpxerr=0
		errdx=0
		dxmdl=twisstwiss.DX[twisstwiss.indx[name]]
		dpxmdl=twisstwiss.DPX[twisstwiss.indx[name]]
		dxp=normal_pro.DPX[normal_pro.indx[name]]
		errdxp=abs(errdmax.DX[errdmax.indx[name]]-errdmin.DX[errdmin.indx[name]])
		
		dpxp=normal_pro.DPX[normal_pro.indx[name]]
		dpxperr="0"
			 
		dxb=back_pro.DX[back_pro.indx[name]]
		errdxb=abs(errdmaxb.DX[errdmaxb.indx[name]]-errdminb.DX[errdminb.indx[name]])
		
		dpxb=-back_pro.DPX[back_pro.indx[name]]
		dpxberr="0"

		datadx.append(name+" "+str(s)+" "+str(dx)+" "+str(errdx)+" "+str(dxmdl)+" "+str(dxp)+" "+str(errdxp)+" "+str(dxb)+" "+str(errdxb))
		datadpx.append(name+" "+str(s)+" "+str(dpx)+" "+str(dpxerr)+" "+str(dpxmdl)+" "+str(dpxp)+" "+str(dpxperr)+" "+str(dpxb)+" "+str(dpxberr))

		# for dy and dpy
		s=twisstwiss.S[twisstwiss.indx[name]]
		dx=0
		dpx=0
		errdx=0
		dpxerr=0
		dxmdl=twisstwiss.DY[twisstwiss.indx[name]]
		dpxmdl=twisstwiss.DPY[twisstwiss.indx[name]]
		dxp=normal_pro.DPY[normal_pro.indx[name]]
		errdxp=abs(errdmax.DY[errdmax.indx[name]]-errdmin.DY[errdmin.indx[name]])

		dpxp=normal_pro.DPY[normal_pro.indx[name]]
		dpxperr="0"
			 
		dxb=back_pro.DY[back_pro.indx[name]]
		errdxb=abs(errdmaxb.DY[errdmaxb.indx[name]]-errdminb.DY[errdminb.indx[name]])

		dpxb=back_pro.DPY[-back_pro.indx[name]]
		dpxberr="0"
		datady.append(name+" "+str(s)+" "+str(dx)+" "+str(errdx)+" "+str(dxmdl)+" "+str(dxp)+" "+str(errdxp)+" "+str(dxb)+" "+str(errdxb))
		datadpy.append(name+" "+str(s)+" "+str(dpx)+" "+str(dpxerr)+" "+str(dpxmdl)+" "+str(dpxp)+" "+str(dpxperr)+" "+str(dpxb)+" "+str(dpxberr))
		

                 #########

	print namename+" is Finished"


print "printing tables to file"
createTables("sbsbetx_"+label+".out",savepath,columnnamesbetx,variablenamebetx,databetx,mainvariablebetx,mainvaluebetx)
createTables("sbsampbetx_"+label+".out",savepath,columnnamesampbetx,variablenameampbetx,dataampbetx,mainvariableampbetx,mainvalueampbetx)
createTables("sbsalfx_"+label+".out",savepath,columnnamesalfx,variablenamealfx,dataalfx,mainvariablealfx,mainvaluealfx)
createTables("sbsbety_"+label+".out",savepath,columnnamesbety,variablenamebety,databety,mainvariablebety,mainvaluebety)
createTables("sbsampbety_"+label+".out",savepath,columnnamesampbety,variablenameampbety,dataampbety,mainvariableampbety,mainvalueampbety)
createTables("sbsalfy_"+label+".out",savepath,columnnamesalfy,variablenamealfy,dataalfy,mainvariablealfy,mainvaluealfy)
createTables("sbsdx_"+label+".out",savepath,columnnames1dx,variablename1dx,datadx,mainvariabledx,mainvaluedx)
createTables("sbsdpx_"+label+".out",savepath,columnnames2dx,variablename2dx,datadpx,mainvariabledx,mainvaluedx)
createTables("sbsdy_"+label+".out",savepath,columnnames1dy,variablename1dy,datady,mainvariabledy,mainvaluedy)
createTables("sbsdpy_"+label+".out",savepath,columnnames2dy,variablename2dy,datadpy,mainvariabledy,mainvaluedy)


print "sbs is finished"




		



	
