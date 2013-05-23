################################################################
#                                                              #
#  @ Glenn Vanbavinckhove  (gvanbavi@cern.ch)=> Date: 11/09/09 #
#                                                              #
################################################################
#
#  !=> SegementBySegment_0.0.py : - Construction of the program (11/09/09)
#  !=> SegementBySegment_0.1.py : - Adding the dispersion to the program (23/09/09)
#                                 - Make distuinsh between BPMs and instuments (24/09/09)
#                                 - Adding output tables (sbs...out) (21/10/09)
#                                 - Adding coupling propogation (21/10/09)
#                                 - changing output when segment is start and end (24/11/09)
#  !=> SegementBySegment_0.22.py :- Fixing some small bugs
#                                 - Added Instruments
#
#  !=> SegementBySegment_0.23.py :- Adding summary list for instruments
#                                 - Applying better cuts for beta
#  !=> SegementBySegment_0.24.py :- Including one/two beam matching (two beam for IP) (7/12/09)
#    
#  !=> SegementBySegment_0.25.py :- Adding phase to the propagation (10/12/09)
#
#  !=> SegementBySegment_0.26.py :- General clean-up (19/01/10)
#                                 - Adding  Total phase to the propagation (19/01/10)
#                                 - Coupling propogation (20/01/10)
#
#  !=> SegementBySegment_0.27.py :- Fixing mistake of total phase (1/03/10)
#                                 - Include real and imaginary part in the output of coupling (1/03/10)
#
#  !=> SegementBySegment_0.28.py : -Coupling initial conditions added to MADX segment  24 March 2010
#                                  -New file getfterms.py required to convert MADX C matrix to observable f terms
#
#



###### imports
from optparse import OptionParser
from metaclass import twiss
import os,sys
from math import sqrt,cos,sin,pi
from datetime import date




###### optionparser
parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                metavar="ACCEL", default="LHCB1",dest="accel")
parser.add_option("-f", "--path", # assumes that output is same as input
                help="Path to measurement files",
                metavar="PATH", default="./", dest="path")
parser.add_option("-i", "--path2", # assumes that output is same as input
                help="Path to second measurement files",
                metavar="PATH2", default="./", dest="path2")
parser.add_option("-s", "--start",
                help="give start,endbpm,name (multiple allowed) eg: start1,end1,name1,start2,end2,name2,...",
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
                help="Output path",
                metavar="SAVE", default="./", dest="SAVE")
parser.add_option("-m", "--mad", # assumes that output is same as input
                help="mad link",
		  metavar="mad", default="", dest="mad")
parser.add_option("-b", "--bbsrouce", # assumes that output is same as input
                help="beta beat source",
                metavar="bb", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/", dest="bb")
parser.add_option("-x", "--take", # assumes that output is same as input
                help="take or create madx 0/1",
                metavar="mad", default="0", dest="madpass")
parser.add_option("-z", "--switch", # assumes that output is same as input
                help="switch one/two beam 1/2",
                metavar="switch", default="1", dest="switch")
parser.add_option("-c", "--cuts", # assumes that output is same as input
                help="cuts for beta,disp,coupling bbxm,bbym,bxe,bye,dxm,dym,dxe,dye,f1001e,f1010e",
                metavar="cuts", default="0.1,0.1,10,10,2,2,0.5,0.5,5,5", dest="cuts")
parser.add_option("-e", "--elementswitch", # assumes that output is same as input
                help="switching between segment(0) or element (1)",
                metavar="holyswitch", default="0", dest="holyswitch")


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


def elementandfilter(beta,ampbeta,disp,phase,couple,startbpm,endbpm,twiss,dispersionswitch,couplingswitch,element):

	##### find bpms between
	bpms=[]
	bpmcut=[]
	bpmstranslate={}
	
	startloc=twiss.S[twiss.indx[startbpm]]
	endloc=twiss.S[twiss.indx[endbpm]]

	all=twiss.NAME
	names=[]
	for name in all:

		if "BPM" in name:
			names.append(name)

	cutswitch=0
		
	if endloc<startloc:

		print "segment cut up in two"
		cutswitch=1

	for name in names:
		loc=twiss.S[twiss.indx[name]]

		if cutswitch==0:

			

			if (loc>startloc and loc<endloc):

				bpms.append(name)
				#print "loc ",startloc,loc,endloc,startbpm,name,endbpm
				bpmstranslate[name]=loc#+twiss.LENGTH


		else:
			
			if (loc>startloc):
				bpms.append(name)
				bpmstranslate[name]=loc
			elif (loc<startloc):
				bpmcut.append(name)
				bpmstranslate[name]=loc+twiss.LENGTH

	if cutswitch==1:
		bpms=[startbpm]+bpms+bpmcut+[endbpm]
		bpmstranslate[startbpm]=startloc
		bpmstranslate[endbpm]=endloc+twiss.LENGTH
	else:
		bpms=[startbpm]+bpms+[endbpm]		
	

	####### structure
	#[x,y,mdlcutx,mdlcuty,errcutx,errcuty]

	bpmsbeta={}
	bpmsb=[]
	bpmsdx={}
	bpmsd=[]
	bpmscouple={}
	bpmsc=[]

	for bpm in bpms:

		try:
                        #beta
			betax=beta[0].BETX[beta[0].indx[bpm]]
			betxmdl=beta[0].BETXMDL[beta[0].indx[bpm]]
			betay=beta[1].BETY[beta[1].indx[bpm]]
			betymdl=beta[1].BETYMDL[beta[1].indx[bpm]]

			alfax=beta[0].ALFX[beta[0].indx[bpm]]
			erralfax=beta[0].ERRALFX[beta[0].indx[bpm]]
			err2alfax=beta[0].STDALFX[beta[0].indx[bpm]]
			alfay=beta[1].ALFY[beta[1].indx[bpm]]
			erralfay=beta[1].ERRALFY[beta[1].indx[bpm]]
			err2alfay=beta[1].STDALFY[beta[0].indx[bpm]]

			erralfax=sqrt(erralfax**2+err2alfax**2)
			erralfay=sqrt(erralfay**2+err2alfay**2)
			
			
			errbetax=beta[0].ERRBETX[beta[0].indx[bpm]]
			stdbetax=beta[0].STDBETX[beta[0].indx[bpm]]
			errbetay=beta[1].ERRBETY[beta[1].indx[bpm]]
			stdbetay=beta[1].STDBETY[beta[1].indx[bpm]]

			errx=sqrt(errbetax**2+stdbetax**2)
			erry=sqrt(errbetay**2+stdbetay**2)

			#ampbeta
			ampbetax=ampbeta[0].BETX[ampbeta[0].indx[bpm]]
			ampbetxmdl=ampbeta[0].BETXMDL[ampbeta[0].indx[bpm]]
			ampbetay=ampbeta[1].BETY[ampbeta[1].indx[bpm]]
			ampbetymdl=ampbeta[1].BETYMDL[ampbeta[1].indx[bpm]]
			
			stdbetax=ampbeta[0].BETXSTD[ampbeta[0].indx[bpm]]
			stdbetay=ampbeta[1].BETYSTD[ampbeta[1].indx[bpm]]
					
	                #phase
			#phasex=phase[0].PHASEX[phase[0].indx[bpm]]
			#phasexmdl=phase[0].PHASEXMDL[phase[0].indx[bpm]]
			#phasey=phase[1].PHASEY[phase[1].indx[bpm]]
			#phaseymdl=phase[1].PHASEYMDL[phase[1].indx[bpm]]
			
			#stdphx=phase[0].STDPHASEX[phase[0].indx[bpm]]
			#stdphy=phase[1].STDPHASEY[phase[1].indx[bpm]]

			#print (betax-betxmdl)/betxmdl,beta[2]

			print errx,beta[4],bpm

			if (betax>0) and (betay>0) and (abs((betax-betxmdl)/betxmdl)<beta[2]) and (abs((betay-betymdl)/betymdl)<beta[3]) and (errx<beta[4]) and (erry<beta[5]) and ((betax-errx)>0) and ((betay-erry)>0) and (errx<betax) and (erry<betay):

				#print "pass"

				bpmsbeta[bpm]=[bpmstranslate[bpm],betax,errx,alfax,erralfax,betay,erry,alfay,erralfay,ampbetax,stdbetax,ampbetay,stdbetay]
				bpmsb.append(bpm)


			else:
				print "BPM ",bpm," didnt pass cuts for beta"

		except:
			print "BPM is not in the datafile for beta",bpm

		try:

			#disperion
			if dispersionswitch==1:
				dispx=disp[0].DX[disp[0].indx[bpm]]
				dispxmdl=disp[0].DXMDL[disp[0].indx[bpm]]
				disppx=disp[0].DPX[disp[0].indx[bpm]]
				disppxmdl=disp[0].DPXMDL[disp[0].indx[bpm]]
				dispy=disp[1].DY[disp[1].indx[bpm]]
				dispymdl=disp[1].DYMDL[disp[1].indx[bpm]]
				dispy=disp[1].DY[disp[1].indx[bpm]]
				dispymdl=disp[1].DYMDL[disp[1].indx[bpm]]

				stddx=ampbeta[0].ERRDX[disp[0].indx[bpm]]
				stddy=ampbeta[1].ERRDY[disp[1].indx[bpm]]
				stddpx=ampbeta[0].ERRDPX[disp[0].indx[bpm]]
				stddpy=ampbeta[1].ERRDPY[disp[1].indx[bpm]]
				

				if(abs(dispx-dispxmdl)<disp[2] and abs(dispy-dispymdl)<disp[3] and  stddx<disp[4] and stddy<disp[5]):
					
					bpmsdx[bpm]=[bpmstranslate[bpm],dispx,stddx,disppx,stddpx,dispy,stddy,disppy,stddpy]
					bpmsd.append(bpm)
				else:

					print "BPM ",bpm," didnt pass cuts for dispersion"

		except:
			print "BPM is not in the datafile for dispersion",bpm

	
		try:
		
	                #coupling
			print "switch",couplingswitch
			#sys.exit()
			if couplingswitch==1:
				#location=couple[0].S[couple[0].indx[bpm]]
				#f1001=couple[0].F1001[couple[0].indx[bpm]]
				#std1001=couple[0].FSTD1[couple[0].indx[bpm]]
				#f1010=couple[0].F1010[couple[0].indx[bpm]]
				#stdf1010=couple[0].FSTD2[couple[0].indx[bpm]]

				location=couple[0].S[couple[0].indx[bpm]]
				f1001w=couple[0].F1001W[couple[0].indx[bpm]]
				std1001=couple[0].FWSTD1[couple[0].indx[bpm]]
				f1001r=couple[0].F1001R[couple[0].indx[bpm]]
				f1001i=couple[0].F1001I[couple[0].indx[bpm]]
				f1010w=couple[0].F1010W[couple[0].indx[bpm]]
				stdf1010=couple[0].FWSTD2[couple[0].indx[bpm]]
				f1010r=couple[0].F1010R[couple[0].indx[bpm]]
				f1010i=couple[0].F1010I[couple[0].indx[bpm]]


				#F1001W FWSTD1 F1001R F1001I F1010W FWSTD2 F1010R F1010I
				print std1001,couple[1],stdf1010,couple[2]
				if(std1001<couple[1] and stdf1010<couple[2]):
					bpmscouple[bpm]=[bpmstranslate[bpm],f1001w,f1010w,std1001,stdf1010,f1001r,f1001i,f1010r,f1010i]
					print "Position "+str(bpmstranslate[bpm])
					#print "adding bpm for coupling"
					bpmsc.append(bpm)
				else:
					print "BPM",bpm," didnt pass cuts for coupling"
		except:
			print "BPM is not in the datafile for coupling",bpm

		#sys.exit()

	if len(bpmsb)<2:

		print "Fatal: Not enough BPM for beta propogation, please increase border => System exit"
		sys.exit()
	else:
		print "The length of the BPMs is ",len(bpmsb)
		

	if (len(bpmsd)<2 and dispersionswitch==1):
		print "Warning: Not enough BPM for dispersion propogation, please increase border. "
		disp=0
	if (len(bpmsc)<2 and couplingswitch==1):
		print "Warning: Not enough BPM for coupling propogation, please increase border"
		couplingswitch=0


	#sys.exit()

	

	if "empty" not in element: # only for element
		locationindx=twiss.indx[element]
		#for beta
		indS=0
		indE=10000000
		bpmsB="empty"
		bpmeB="empty"
		for bpm in bpmsb:
			#print bpm
			indx=twiss.indx[bpm]
		
			if locationindx>indx and indx>indS:
				indxS=indx
				bpmsB=bpm			
			if locationindx<indx and indx<indE:
				indxE=indx
				bpmeB=bpm
				break

	
		indS=0
		indE=10000000
		bpmsD="empty"
		bpmeD="empty"
		for bpm in bpmsd:
			indx=twiss.indx[bpm]
			if locationindx>indx and indx>indS:
				indxS=indx
				bpmsD=bpm
			if locationindx<indx and indx<indE:
				indxE=indx
				bpmeD=bpm
				break

		indS=0
		indE=10000000
		bpmsC="empty"
		bpmeC="empty"
		for bpm in bpmsc:
			indx=twiss.indx[bpm]
			if locationindx>indx and indx>indS:
				indxS=indx
				bpmsC=bpm
			if locationindx<indx and indx<indE:
				indxE=indx
				bpmeC=bpm
				break

		

		elbpms=[bpmsB,bpmeB,bpmsD,bpmeD,bpmsC,bpmeC]

		print bpmsB,bpmeB

	else:
		elbpms=["empty","empty","empty","empty","empty","empty"]

	#sys.exit()
	
	return [bpmsb,bpmsbeta],[bpmsd,bpmsdx],[bpmsc,bpmscouple],disp,couplingswitch,elbpms


def getTwiss(filee,element):

    print "loading file ",filee

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

def run4mad(path,hor,ver,hore,vere,dp,dpe,startbpm,endbpm,name, fs, exppath):


    ##
    #  Copy getfterms.py locally to be used by MADx
    cpath=options.bb
    os.system('cp '+cpath+'/SegmentBySegment/getfterms.py'+' '+path+'/')
    
    if options.accel=="LHCB2":

        dire=-1
        start="MKI.A5R8.B2"
        beam="B2"
        
    elif options.accel=="LHCB1":

        dire=1
	start="MSIA.EXIT.B1"   #  compatible with repository
#	start="IP2"   #  compatible with repository
        beam="B1"

   

    ### check on error propogation
    errbetx=hor[1]
    betx=hor[0]
    #if (hor[0]-hor[1])<0:errbetx=0;betx=hor[0]
    errbety=ver[1]
    bety=ver[0]
    #if (ver[0]-ver[1])<0:errbety=0;bety=ver[0]
    errbetxb=hore[1]
    betxb=hore[0]
    #if (hore[0]-hore[1])<0:errbetxb=0;betxb=hore[0]
    errbetyb=vere[1]
    betyb=vere[0]    
    #if (vere[0]-vere[1])<0:errbetyb=0;betyb=0
    f1001r=fs[0]
    f1001i=fs[1]
    f1010r=fs[2]
    f1010i=fs[3]

    
     
    filename=path+'/var4mad.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%BETX/\''+str(betx)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETY/\''+str(bety)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETX/\''+str(errbetx)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRBETY/\''+str(errbety)+'\'/g\' \\\n')
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
    
    file4nad.write('    -e \'s/%ENDBX/\''+str(betxb)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDBY/\''+str(betyb)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRENDBX/\''+str(errbetxb)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ERRENDBY/\''+str(errbetyb)+'\'/g\' \\\n')
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
    
    file4nad.write('    -e \'s/%STARTFROM/\''+str(startbpm)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDAT/\''+str(endbpm)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(name)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DIRE/\''+str(dire)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%START/\''+str(start)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BEAM/\''+str(beam)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PATH/\'\"'+str(path.replace('/','\/'))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001R/\''+str(f1001r)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1001I/\''+str(f1001i)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1010R/\''+str(f1010r)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%F1010I/\''+str(f1010i)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%EXP/\''+exppath.replace("/","")+'\'/g\' \\\n')
    
    
    
    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/job.InterpolateBetas.testCoupl.mask > '+path+'/t_'+str(name)+'.madx \n')

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))

    runmad(path,name)

   
    
def runmad(path,name):
	
	os.system(options.mad+'madx < '+path+'t_'+str(name)+'.madx')
	
   
  
def run4plot(path,spos,epos,beta4plot,cpath,meapath,name):

    filename=path+'/var4plot.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%PATH/\'\"'+str(path.replace("/","\/"))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%EndPoint/\''+str(epos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%StartPoint/\''+str(spos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(name)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETA/\''+str(beta4plot)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%MEA/\'\"'+str(meapath.replace("/","\/"))+'\"\'/g\' \\\n')    
    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.coupl.mask > '+path+'/gplot_'+name)

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))
    #os.system("chmod 777 "+str(path+'/gplot_'+name))
    os.system("gnuplot "+str(path+'/gplot_'+name))

 
   
def GetPhaseEM(exp, mod):
    phasem=[]
    bpm1=[]
    bpm2=[]
    s1=[]
    s2=[]
    phaseexp=[]
    all=mod.NAME
    names=[]

    for name in all:

        if "BPM" in name:

	    names.append(name)

    
    for elm in names:
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

def reversetable(path,name):
	file=open(path+"/twiss_"+name+"_back_rev.dat",'w')
	base=twiss(path+"/twiss_"+name+"_back.dat")

	file.write("* NAME                                S               BETX               ALFX               BETY               ALFY    DX       DPX     DY     DPY   MUX   MUY\n")
	file.write("$ %s                                %le                %le                %le                %le                %le      %le                %le      %le   %le            %le                %le \n")
	
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
	mux=base.MUX
	muy=base.MUY

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
		muxx=mux[index]
		muyy=muy[index]

		#print bpm, dex

		file.write(str(bpm)+' '+str(endpos-ss)+' '+str(bex)+' '+str(alx)+' '+str(bey)+' '+str(aly)+'  '+str(dex)+' '+str(depx)+' '+str(dey)+' '+str(depy)+' '+str(muxx)+' '+str(muyy)+'\n')
		
	file.close()


def rotateparts(datacouple,twiss):

	bpms=datacouple[0]
	data=datacouple[1]
	startbpm=bpms[0]
	endbpm=bpms[len(bpms)-1]

	twiss.Cmatrix()

	
	coupleddata={}
	xx=1j

	#[data[bpms[bpm+1]][0],f1001real,f1001imag,f1001stdr,f1001basic,stdf1001basic,f1001basicb,stdf1001basicb,f1001M,f1010real,f1010imag,f1010stdr,f1010M,f1010basic,stdf1010basic,f1010basicb,stdf1010basicb]
	
	f1001basic=data[startbpm][5]+data[startbpm][6]*xx
	f1010basic=data[startbpm][7]+data[startbpm][8]*xx
	stdf1001basic=data[startbpm][3]
	stdf1010basic=data[startbpm][4]


	f1001basicb=data[endbpm][5]+data[endbpm][6]*xx
	f1010basicb=data[endbpm][7]+data[endbpm][8]*xx
	stdf1001basicb=data[endbpm][3]
	stdf1010basicb=data[endbpm][4]


	coupleddata[startbpm]=[data[startbpm][0],data[startbpm][5],data[startbpm][6],data[startbpm][3],f1001basic,data[startbpm][3],f1001basic,data[startbpm][3],twiss.F1001W[twiss.indx[startbpm]],data[startbpm][7],data[startbpm][8],data[startbpm][4],twiss.F1010W[twiss.indx[startbpm]],f1010basic,data[startbpm][4],f1010basic,data[startbpm][4]]
	coupleddata[endbpm]=[data[endbpm][0],data[endbpm][5],data[endbpm][6],data[endbpm][3],f1001basicb,data[endbpm][3],f1001basicb,data[endbpm][3],twiss.F1001W[twiss.indx[endbpm]],data[endbpm][7],data[endbpm][8],data[endbpm][4],twiss.F1010W[twiss.indx[endbpm]],f1010basicb,data[endbpm][4],f1010basicb,data[endbpm][4]]
	
	for bpm in range(len(bpms)-1):

		#print bpms[bpm+1]
		#sys.exit()

		f1001r=data[bpms[bpm+1]][1]
		f1001stdr=data[bpms[bpm+1]][3]
		f1010r=data[bpms[bpm+1]][2]
		f1010stdr=data[bpms[bpm+1]][4]

		
		f1001real=data[bpms[bpm+1]][5]
		f1001imag=data[bpms[bpm+1]][6]
		f1010real=data[bpms[bpm+1]][7]
		f1010imag=data[bpms[bpm+1]][8]
		


		phx12=twiss.MUX[twiss.indx[bpms[bpm+1]]]-twiss.MUX[twiss.indx[bpms[bpm]]]
		phy12=twiss.MUY[twiss.indx[bpms[bpm+1]]]-twiss.MUY[twiss.indx[bpms[bpm]]]

		delta1001=(phx12-phy12)*2*pi
		delta1010=(phx12+phy12)*2*pi
		#print len(bpms)-1-bpm,bpms[len(bpms)-1-bpm]
	        phx12=-twiss.MUX[twiss.indx[bpms[bpm]]]+twiss.MUX[twiss.indx[bpms[bpm+1]]]
		phy12=-twiss.MUY[twiss.indx[bpms[bpm]]]+twiss.MUY[twiss.indx[bpms[bpm+1]]]

		delta1001b=-(phx12-phy12)*2*pi
		delta1010b=-(phx12+phy12)*2*pi

		#print cos(delta1001),sin(delta1001),abs(f1001basic),f1001basic
		
		f1001basic=complex(f1001basic.real,f1001basic.imag)*complex(cos(delta1001),sin(delta1001))
		f1010basic=complex(f1010basic.real,f1010basic.imag)*complex(cos(delta1001),sin(delta1001))
		stdf1001basic=stdf1001basic*complex(cos(delta1010),sin(delta1010))	       
		stdf1010basic=stdf1010basic*complex(cos(delta1010),sin(delta1010))

		f1001basicb=complex(f1001basicb.real,f1001basicb.imag)*complex(cos(delta1001b),sin(delta1001b))
		f1010basicb=complex(f1010basicb.real,f1010basicb.imag)*complex(cos(delta1001b),sin(delta1001b))
		stdf1001basicb=stdf1001basicb*complex(cos(delta1010b),sin(delta1010b))	       
		stdf1010basicb=stdf1010basicb*complex(cos(delta1010b),sin(delta1010b))

		f1001M=twiss.F1001W[twiss.indx[bpms[bpm+1]]]
		f1010M=twiss.F1010W[twiss.indx[bpms[bpm+1]]]

		print f1001real,f1001basic.real,f1001basicb.real,cos(delta1001),sin(delta1001)

		
		coupleddata[bpms[bpm+1]]=[data[bpms[bpm+1]][0],f1001real,f1001imag,f1001stdr,f1001basic,abs(stdf1001basic),f1001basicb,abs(stdf1001basicb),f1001M,f1010real,f1010imag,f1010stdr,f1010M,f1010basic,abs(stdf1010basic),f1010basicb,abs(stdf1010basicb)]

	#sys.exit()
	
	return [bpms,coupleddata]


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

	print len(columnnames),len(paranames),outputname
	

	if len(columnnames)==len(paranames):


		filefile.write('* '+'  '.join(columnnames)+'\n')
		filefile.write('$ '+'  '.join(paranames)+'\n')

		#writing data
		#
		# data[x][y]
		#
		#   x = data set
		#   y = data in data set
		#   => has to bed added in same order
		#
		
		for y in range(len(data)):
		
			filefile.write(str(data[y])+'\n')


	
	
	else:

		print "cannot create table names for columns are not equal => system exit "
		sys.exit()

	filefile.close()



def filterbpms(filex,filey,bpmstart,bpmend,listt,listcom,accel):

     bpms=[]
     names=filex.NAME
     go=0

     for name in names:

	     try: # if both planes contain the name

		  
	             bx=filex.BETX[filex.indx[name]]
	             by=filey.BETY[filey.indx[name]]
		     errbetx=(sqrt(filex.ERRBETX[filex.indx[name]]**2+filex.STDBETX[filex.indx[name]]**2)/filex.BETX[filex.indx[name]])*100	     
		     errbety=(sqrt(filey.ERRBETY[filey.indx[name]]**2+filey.STDBETY[filey.indx[name]]**2)/filey.BETY[filey.indx[name]])*100
		     print bpmstart,name
		     if bpmstart in name:
			     go=1
			     print "found start"
			     #sys.exit()
		     if go==1:
			     if ((bx and by) >0) and ((errbetx and errbety) <20):

				     bpms.append(name)

			     else:
				     print "not excepting bpm because of high error : "+name
		     if bpmend==name:
			     break

	     except:
	     
		     print "Beta not found for both planes : "+name
     #print len(bpms)
     #print 	bpms[0],bpms[len(bpms)-1]	   
     #sys.exit()		     
     correctors=correctorselection(bpms[0],bpms[len(bpms)-1],listt,listcom,accel)

     return bpms,correctors

def correctorselection(Bpmstart,Bpmend,listt,listcom,accel):
   
    names=twiss(listt).NAME
    if accel=="LHCB1":
	    names=names+twiss(listcom).NAME
    go=0
    correctorlist=[]
    cleanerl=Bpmstart.split(".")
    cleanerr=Bpmend.split(".")
    cleanerl2=list(cleanerl[1])
    onefilter=cleanerl2[len(cleanerl2)-2]+cleanerl2[len(cleanerl2)-1]
    cleanerr2=list(cleanerr[1])
    twofilter=cleanerr2[len(cleanerr2)-2]+cleanerr2[len(cleanerr2)-1]
    for name in names:
        #print name.upper()
	namess=name.split('.')
	#print namess
	if onefilter in namess[1].upper() or twofilter in namess[1].upper():
	        #print name
		correctorlist.append(name)

    #sys.exit()
    return correctorlist
			





def creatematching(betaxb1,betayb1,betaxb2,betayb2,path,beamswitch,bpmstart,bpmend,listt1,listt2,listcom,name,accel1,accel2):

    matcher=open(path+'/matchjob4_'+name+'_sbs.madx','w')
    #####
    #	    # start        #
    #####
    if accel1=="LHCB1":
	    bpmstart2=bpmstart.replace("B1","B2")
	    bpmend2=bpmend.replace("B1","B2")
	    accel2="LHCB2"
	    accel1="LHCB1"
    else:
	    bpmstart2=bpmstart.replace("B2","B1")
	    bpmend2=bpmend.replace("B2","B1")
	    accel2="LHCB1"
	    accel1="LHCB2"
    if beamswitch=='2':
	    matcher.write('seqedit, sequence='+accel2+';\n')
	    matcher.write('flatten;\n')
	    matcher.write('cycle, start=MKI.A5R8.B2;;\n')   ###### starting matching
	    matcher.write('endedit;\n')
		    
	    matcher.write('use, period='+accel2+';\n')
	    matcher.write('seqedit, sequence='+accel2+';\n')
	    matcher.write('flatten;\n')
	    matcher.write('cycle, start=MKI.A5R8.B2;\n')   ###### starting matching
	    matcher.write('endedit;\n')
	    matcher.write('use, period='+accel2+', range='+bpmstart2+'/'+ bpmend2+';\n')
    try:
	    bpms1,correctorsb1=filterbpms(betaxb1,betayb1,bpmstart,bpmend,listt1,listcom,accel1)

    except:
	    print "cannot find bpms YET for segment which is cut in two"
	    bpms1=[]
	    correctorsb1=[]
    matcher.write('assign, echo="'+path+'/matching_before.dat";\n')

	    
    for corrector in correctorsb1:
	    matcher.write('value,'+corrector +';\n')
    if beamswitch=='2':
	    try:
		    bpms2,correctorsb2=filterbpms(betaxb2,betayb2,bpmstart2,bpmend2,listt2,listcom,accel2)
	    except:
		    print "cannot find bpms YET for segment which is cut in two"
		    bpms2=[]
		    correctorsb2=[]
	    for corrector in correctorsb2:
		    matcher.write('value,'+corrector +';\n')
    matcher.write('assign, echo=terminal;\n')
    matcher.write('match, use_macro;\n')
	    
    #####
    #
    # constraint
    #
    #####
    if len(bpms1)>0:
	    name1=bpms1[0]
	    bx1=betaxb1.BETX[betaxb1.indx[name1]]
	    ax1=betaxb1.ALFX[betaxb1.indx[name1]]
	    by1=betayb1.BETY[betayb1.indx[name1]]
	    ay1=betayb1.ALFY[betayb1.indx[name1]]
	    bpms1.remove(name1) # not placing the first bpm in the constraint
    ##### variation list beam 1
    for corrector in correctorsb1:
	    matcher.write('vary, name='+corrector+';\n')	   
    if beamswitch=='2':
	    if len(bpms1)>0:
		    name2=bpms2[0]
		    bx2=betaxb2.BETX[betaxb2.indx[name2]]
		    ax2=betaxb2.ALFX[betaxb2.indx[name2]]
		    by2=betayb2.BETY[betayb2.indx[name2]]
		    ay2=betayb2.ALFY[betayb2.indx[name2]]
		    bpms2.remove(name2)
            ##### variation list beam 2
	    for corrector in correctorsb2:
		    matcher.write('vary, name='+corrector+';\n')
    if len(bpms1)>0:
	    matcher.write('M1: MACRO={ \n')
	    matcher.write('select, flag=twiss,pattern="bpm", column=name, s, betx, bety;\n')
	    matcher.write('twiss, sequence='+accel1+', betx='+str(bx1)+', alfx='+str(ax1)+', bety='+str(by1)+', alfy='+str(ay1)+',dx=0,dy=0,dpx=0,dpy=0, table=t1;\n') # beam1
	    
	    if beamswitch=='2':
		    matcher.write('twiss, sequence='+accel2+', betx='+str(bx2)+', alfx='+str(ax2)+', bety='+str(by2)+', alfy='+str(ay2)+',dx=0,dy=0,dpx=0,dpy=0, table=t2;\n') # beam2
	    matcher.write('}\n')
    #####
    #
    # beam1
    #
    #####
 
    for name in bpms1:

	    bx=betaxb1.BETX[betaxb1.indx[name]]
	    by=betayb1.BETY[betayb1.indx[name]]
		
	    matcher.write('CONSTRAINT, EXPR= table(t1,'+name+',bety ) ='+str(by)+' ;\n')
	    matcher.write('CONSTRAINT, EXPR= table(t1,'+name+',betx ) ='+str(bx)+' ;\n')


    #####
    #
    # beam2
    #
    #####
    if beamswitch=='2':
	    for name in bpms2:

		    bx=betaxb2.BETX[betaxb2.indx[name]]
		    by=betayb2.BETY[betayb2.indx[name]]

		    matcher.write('CONSTRAINT, EXPR= table(t2,'+name+',bety ) ='+str(by)+' ;\n')
		    matcher.write('CONSTRAINT, EXPR= table(t2,'+name+',betx ) ='+str(bx)+' ;\n')
    #####
    #
    # ending
    #
    #####
    matcher.write('simplex, tolerance=1e-12;\n'+
	                  'endmatch;\n')    
    #####
    #
    # writing to file
    #
    #####
    matcher.write('assign, echo="'+path+'/matching_output.dat";\n')
    for corrector in correctorsb1:
	    matcher.write('value,'+corrector +';\n')
    if beamswitch=='2':
	    for corrector in correctorsb2:
		    matcher.write('value,'+corrector +';\n')
    matcher.write('assign, echo=terminal;\n')
    
    matcher.close()


####### 
#
# #
# main part #
# #
#
######
#file4seg=twiss(options.bb+"/MODEL/"+options.accel+"/sbslist/"+options.segf)
path=options.path
beamswitch=options.switch
#if options.accel=="LHCB1":
#	accel2="LHCB2"
#	accel1="LHCB1"
#	listt2=options.bb+'/MODEL/'+accel2+'/sbslist//listb2.dat'
#	listt1=options.bb+'/MODEL/'+accel1+'/sbslist/listb1.dat'
#else:
#	accel2="LHCB1"
#	accel1="LHCB2"
#	listt2=options.bb+'/MODEL/'+accel1+'/sbslist//listb2.dat'
#	listt1=options.bb+'/MODEL/'+accel2+'/sbslist/listb1.dat'

#listcom=options.bb+'/MODEL/'+accel1+'/sbslist/listcom.dat'
##########
#
# LOADING DATA
#
##########

#
# beam1
#
filedatax=twiss(path+"/getbetax.out")
betxtwiss=filedatax
filedatay=twiss(path+"/getbetay.out")
betytwiss=filedatay
filedataxA=twiss(path+"/getampbetax.out")
ampbetxtwiss=filedataxA
filedatayA=twiss(path+"/getampbetay.out")
ampbetytwiss=filedatayA
filephasex=twiss(path+"/getphasey.out")
filephasey=twiss(path+"/getphasex.out")
filephasextot=twiss(path+"/getphasetotx.out")
filephaseytot=twiss(path+"/getphasetoty.out")

	
if beamswitch=='2':
	try:
		filedatax2=twiss(options.path2+"/getbetax.out")
		filedatay2=twiss(options.path2+"/getbetay.out")
		print "second beam found"
		#sys.exit()
	except:
		print "you selected two beam match ... but for path2 the files doesn't exist"
		beamswitch='1'

		#sys.exit()
		
### check if dispersion exist
try:
	filedx=twiss(path+"/getDx.out")
	filedy=twiss(path+"/getDy.out")
	disp=1
except:
	print "no dispersion files... will continue without taking into account dispersion"
	filedx=[]
	filedy=[]
	disp=0

### check if coupling exist
try:
	filecouple=twiss(path+"/getcouple.out")
	coupleswitch=1
except:
	print "no coupling file... will continue without taking into account coupling"
	filecouple=[]
	coupleswitch=0
	fs=[0,0,0,0]



#
# beam2
#
if beamswitch=='2':
	try:
		filedatax2=twiss(options.path2+"/getbetax.out")
		filedatay2=twiss(options.path2+"/getbetay.out")
		print "second beam found"
	except:
		print "you selected two beam match ... but for path2 the files doesn't exist"
		filedatax2=twiss(options.path+"/getbetax.out")
		filedatay2=twiss(options.path+"/getbetay.out")
		beamswitch='1'
else: # silly trick
		filedatax2=twiss(options.path+"/getbetax.out")
		filedatay2=twiss(options.path+"/getbetay.out")	


# checking if end element is BPM or instrument (IP,collimators)
elementswitch=-1
savepath=options.SAVE
list2run=options.segf.split(',')
cuts=options.cuts.split(',')
start={}
end={}
names=[]

twisspath=options.twiss
if twisspath=="./":
	twisspath=options.bb+"/MODEL/"+options.accel+"/nominal.opt/twiss.dat"
	
twisstwiss=twiss(twisspath)

	
if options.holyswitch=="0":
	if (len(list2run)/3)!=1:
		print "You didnt gave the right input for segment start,end,name => System exit"
		sys.exit()
	for run in range(len(list2run)-2):
		
		start[list2run[run+2]]=list2run[run]
		end[list2run[run+2]]=list2run[run+1]
		names.append(list2run[run+2])


else:
	for run in range(len(list2run)):
		start[list2run[run]]=twisstwiss.NAME[0]
		end[list2run[run]]=twisstwiss.NAME[len(twisstwiss.NAME)-1]
		names.append(list2run[run])


####### for data files
mainvariablebetx=['NFILES']
mainvaluebetx=[filedatax.COUNT[0]]
columnnamesbetx=['NAME','S','BETX','ERRBETX','BETXMDL','BETXP','ERRBETXP','BETXB','ERRBETXB']
variablenamebetx=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
databetx=[]
# for ampbetx
mainvariableampbetx=['NFILES']
mainvalueampbetx=[filedatax.COUNT[0]]
columnnamesampbetx=['NAME','S','BETX','ERRBETX','BETXMDL','BETXP','ERRBETXP','BETXB','ERRBETXB']
variablenameampbetx=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
dataampbetx=[]
# for alfx
columnnamesalfx=['NAME','S','ALFX','ERRALFX','ALFXMDL','ALFXP','ERRALFXP','ALFXB','ERRALFXB']
variablenamealfx=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablealfx=['NFILES']
mainvaluealfx=[filedatax.COUNT[0]]
dataalfx=[]
# for bety
betytwiss=twiss(path+"/getbetay.out")
columnnamesbety=['NAME','S','BETY','ERRBETY','BETYMDL','BETYP','ERRBETYP','BETYB','ERRBETYB']
variablenamebety=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablebety=['NFILES']
mainvaluebety=[filedatay.COUNT[0]]
databety=[]
# for ampbety
columnnamesampbety=['NAME','S','BETY','ERRBETY','BETYMDL','BETYP','ERRBETYP','BETYB','ERRBETYB']
variablenameampbety=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariableampbety=['NFILES']
mainvalueampbety=[filedatay.COUNT[0]]
dataampbety=[]
# for alfy
columnnamesalfy=['NAME','S','ALFY','ERRALFY','ALFYMDL','ALFYP','ERRALFYP','ALFYB','ERRALFYB']
variablenamealfy=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablealfy=['NFILES']
mainvaluealfy=[filedatay.COUNT[0]]
dataalfy=[]
#for dx and dpx
dx4twiss=path+"/getDx.out"
columnnames1dx=['NAME','S','DX','ERRDX','DXMDL','DXP','ERRDXP','DXB','ERRDXB']
variablename1dx=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
columnnames2dx=['NAME','S','DPX','ERRDPX','DPXMDL','DPXP','ERRDPXP','DPXB','ERRDPXB']
variablename2dx=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariabledx=['NFILES']
mainvaluedx=[0]
datadx=[]
datadpx=[]
# for dy and dpy
dy4twiss=path+"/getDy.out"
columnnames1dy=['NAME','S','DY','ERRDY','DYMDL','DYP','ERRDYP','DYB','ERRDYB']
variablename1dy=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
columnnames2dy=['NAME','S','DPY','ERRDPY','DPYMDL','DPYP','ERRDPYP','DPYB','ERRDPYB']
variablename2dy=['%s','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariabledy=['NFILES']
mainvaluedy=[0]
datady=[]
datadpy=[]
#for phasextot
columnnamesphasextot=['NAME','S','PHASEX','ERRPHASEX','PHASEXP','ERRPHASEXP','MdiffP','PHASEXB','ERRPHASEXB','MdiffB']
variablenamephasextot=['%s','%le','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablephasextot=['NFILES']
mainvaluephasextot=[filephasextot.COUNT[0]]
dataphasextot=[]
#for phaseytot
columnnamesphaseytot=['NAME','S','PHASEY','ERRPHASEY','PHASEYP','ERRPHASEYP','MdiffP','PHASEYB','ERRPHASEXB','MdiffB']
variablenamephaseytot=['%s','%le','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablephaseytot=['NFILES']
mainvaluephaseytot=[filephaseytot.COUNT[0]]
dataphaseytot=[]
# for couple
columnnamescouple=['NAME','S','F1001','STDF1001','ReF1001','ImF1001','F1001P','ReF1001P','ImF1001P','STDF1001P','F1001M','F1010','STDF1010','ReF1010','ImF1010','F1010P','ReF1010P','ImF1010P','STDF1010P','F1010M','F1001B','ReF1001B','ImF1001B','STDF1001B','F1010B','ReF1010B','ImF1010B','STDF1010B']
variablenamecouple=['%s','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le','%le']
mainvariablecouple=[]
mainvaluecouple=[]
coupledata=[]


######## big loop

####### Difference between segment or instrument #######
# For segment => Startbpm,endbpm,NameOfSegment
# For instrument => StartBPM,endBpm, instrument (must be excisting name)
####### 

for namename in names:

	#getting common bpms in model and measurment
	#bpms=intersect([filedatax,filedatay])
	#commonbpms=modelIntersect(bpms,twisstwiss)

	passs=1
 	
	if options.holyswitch=="0":
		elementswitch=0
		print "Segment has been choosen"
		databeta,datadx,datacouple,disp,coupleswitch,elbpms=elementandfilter(
			[betxtwiss,betytwiss,cuts[0],cuts[1],cuts[2],cuts[3]],
			[ampbetxtwiss,ampbetytwiss],[filephasextot,filephaseytot],
			[filedx,filedy,cuts[4],cuts[5],cuts[6],cuts[7]],
			[filecouple,cuts[8],cuts[9]],start[namename],
			end[namename],
			twisstwiss,
			disp,coupleswitch,"empty")

	else:
		elementswitch=1
		print "You have choosen an instrument"
		databeta,datadx,datacouple,disp,coupleswitch,elbpms=elementandfilter([betxtwiss,betytwiss,cuts[0],cuts[1],cuts[2],cuts[3]],[ampbetxtwiss,ampbetytwiss],[filephasextot,filephaseytot],[filedx,filedy,cuts[4],cuts[5],cuts[6],cuts[7]],[filecouple,cuts[8],cuts[9]],start[namename],end[namename],twisstwiss,disp,coupleswitch,namename)



	if elementswitch==0:

		startbpm=databeta[0][0]
		endbpm=databeta[0][len(databeta[0])-1]
	else:
		startbpm=elbpms[0]
		endbpm=elbpms[1]
	


	if ("empty" not in startbpm and "empty" not in endbpm):
		hor=[databeta[1][startbpm][1],databeta[1][startbpm][2],databeta[1][startbpm][3],databeta[1][startbpm][4]]
		ver=[databeta[1][startbpm][5],databeta[1][startbpm][6],databeta[1][startbpm][7],databeta[1][startbpm][8]]
		hore=[databeta[1][endbpm][1],databeta[1][endbpm][2],databeta[1][endbpm][3],databeta[1][endbpm][4]]	
		vere=[databeta[1][endbpm][5],databeta[1][endbpm][6],databeta[1][endbpm][7],databeta[1][endbpm][8]]
	else:
		print "WARN : haven't find the appropaites boundaries skipping ",namename
		passs=0
	if disp==1:

		if elementswitch==0:

			startbpm=datadx[0][0]
			endbpm=datadx[0][len(databeta[0])-1]
		else:
			startbpm=elbpms[2]
			endbpm=elbpms[3]

		if ("empty" not in startbpm and "empty" not in endbpm):
		
			dp=[datadx[1][startbpm][1],datadx[1][startbpm][3],datadx[1][startbpm][5],datadx[1][startbpm][7],datadx[1][startbpm][2],datadx[1][startbpm][6]]
			dpe=[datadx[1][endbpm][1],datadx[1][endbpm][3],datadx[1][endbpm][5],datadx[1][endbpm][7],datadx[1][endbpm][2],datadx[1][endbpm][6]]
		else:
			print "WARN : haven't find the appropaites boundaries skipping ",namename
			disp==0
			dp=[0,0,0,0,0,0]
			dpe=[0,0,0,0,0,0]
			
	else:
		dp=[0,0,0,0,0,0]
		dpe=[0,0,0,0,0,0]
		print "No dispersion"

	if coupleswitch==1:
		coupling=rotateparts(datacouple,twisstwiss)
		print len(coupling[1][coupling[0][0]]),len(coupling[0])
		data=datacouple[1]
		f1001r=data[startbpm][5]
		f1001i=data[startbpm][6]
		f1010r=data[startbpm][7]
		f1010i=data[startbpm][8]
		fs=[f1001r,f1001i,f1010r,f1010i]
		#sys.exit()
	else:
		coupling=[[],[]]
		fs=[0,0,0,0]
		print "No coupling"
		
	
	
        #creatematching(filedatax,filedatay,filedatax2,filedatay2,path,beamswitch,element,eelement,listt1,listt2,listcom,label,accel1,accel2)

	if passs==1:
		
		if str(options.madpass)=="0":
			run4mad(savepath,hor,ver,hore,vere,dp,dpe,startbpm,endbpm,namename, fs, options.path)

		else:
			runmad(savepath,namename)
			print "Just rerunning mad"
		reversetable(savepath,namename)


	##################
	#=> switch only for element - to - element
	if elementswitch==0:
	##################
		[hor,ver,dp]=getTwiss(savepath+"/StartPoint.twiss",startbpm)
		startpos=hor[2]
		[hor,ver,dp]=getTwiss(savepath+"/twiss_"+str(namename)+".dat",startbpm)
		betx=hor[0]
		bety=ver[0]
		ax=hor[1]
		ay=ver[1]
		dx=dp[0]
		dpx=dp[1]
		dy=dp[2]
		dpy=dp[3]
		[hor,ver,dp]=getTwiss(savepath+"/StartPoint.twiss",endbpm)
		endpos=hor[2]

		[hor,ver,dp]=getTwiss(savepath+"/twiss.b+.dat",startbpm)
		bplusx=hor[0]
		bplusy=ver[0]
		[hor,ver,dp]=getTwiss(savepath+"/twiss.b-.dat",startbpm)
		bminx=hor[0]
		bminy=ver[0]

		[hor,ver,dp]=getTwiss(savepath+"/twiss.a+.dat",startbpm)
		aplusx=hor[1]
		aplusy=ver[1]
		[hor,ver,dp]=getTwiss(savepath+"/twiss.a-.dat",startbpm)
		aminx=hor[1]
		aminy=ver[1]

		berrx=bplusx-bminx
	        aerrx=aplusx-aminx
		berry=bplusy-bminy
		aerry=aplusy-aminy

		phx=twiss(path+'/getphasex.out')
		phy=twiss(path+'/getphasey.out')

		m=twiss(savepath+"/twiss_"+namename+".dat")

		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phx, m)
		writePhase(savepath+"/phasexEM.out",bpm1, bpm2, s1, s2, phaseexp, phasem)
		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phy, m)
		writePhase(savepath+"/phaseyEM.out",bpm1, bpm2, s1, s2, phaseexp, phasem)

		m=twiss(savepath+"/twiss_"+namename+"_play.dat")

		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phx, m)
		writePhase(savepath+"/phasexEM_play.out",bpm1, bpm2, s1, s2, phaseexp, phasem)
		bpm1, bpm2, s1, s2, phaseexp, phasem = GetPhaseEM(phy, m)
		writePhase(savepath+"/phaseyEM_play.out",bpm1, bpm2, s1, s2, phaseexp, phasem)


		#m=options.twiss
		[hor,ver,dp]=getTwiss(twisspath,startbpm)
		beta4plot=hor[2]

	        ######################
		if options.gra=='1':
			run4plot(savepath,startpos,endpos,beta4plot,options.bb,path,namename)


		################
	        #fileresul=open(savepath+'/resul_'+namename+'.tfs','w')
	        #fileresul.write('* NAME   S     BETX    ERRBETX    ALFX   ERRALFX   BETY   ERRBETY  ALFY   ERRALFX  DX  DPX   DY    DPY\n')
	        #fileresul.write('* %s     %le    %le   %le         %le    %le       %le    %le      %le    %le     %le   %le    %le    %le\n')
	        #fileresul.write(str(endbpm)+' '+str(endpos)+' '+str(betx)+' '+str(berrx)+' '+str(ax)+' '+str(aerrx)+' '+str(bety)+' '+str(berry)+' '+str(ay)+' '+str(aerry)+' '+str(dx)+' '+str(dpx)+' '+str(dy)+' '+str(dpy))
	        #fileresul.close()


	        #### creating tables
		#
		# => opening all the results
		#


		# results from propagation
		normal_pro=twiss(savepath+'/twiss_'+namename+'.dat')
		back_pro=twiss(savepath+'/twiss_'+namename+'_back_rev.dat')
		bpmsbetx=intersect([filedatax,normal_pro,back_pro])
		ampbpmsbetx=intersect([ampbetxtwiss,normal_pro,back_pro])
		bpmsbety=intersect([filedatay,normal_pro,back_pro])
		ampbpmsbety=intersect([ampbetytwiss,normal_pro,back_pro])


		bpmsphasex=intersect([filephasextot,normal_pro,back_pro])
		bpmsphasey=intersect([filephaseytot,normal_pro,back_pro])
			
		if disp==1:

			Dx=twiss(path+"/getDx.out")
			bpmsdx=intersect([Dx,normal_pro,back_pro])
			Dy=twiss(path+"/getDy.out")
			bpmsdy=intersect([Dy,normal_pro,back_pro])
			try:
				mainvaluedy=[Dy.COUNT[0]]
				mainvaluedx=[Dx.COUNT[0]]
			except:
				mainvaluedy=[0]
				mainvaluedx=[0]
	
		else:
			print "no dispersion file"
			bpmsdx=[]
			mainvaluedx=[0]
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
		
		print "filling table for ampbetx"
		for bpm in ampbpmsbetx:
			name=bpm[1]
			s=ampbetxtwiss.S[ampbetxtwiss.indx[name]]
			#if ampbetxtwiss.S[ampbetxtwiss.indx[name]]<ampbetxtwiss.S[ampbetxtwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=ampbetxtwiss.BETX[ampbetxtwiss.indx[name]]
			errbet=ampbetxtwiss.BETXSTD[ampbetxtwiss.indx[name]]
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
			#if betxtwiss.S[betxtwiss.indx[name]]<betxtwiss.S[betxtwiss.indx[element]]:s=s+twisstwiss.LENGTH
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
			#if betytwiss.S[betytwiss.indx[name]]<betytwiss.S[betytwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=betytwiss.BETY[betytwiss.indx[name]]
			errbet=sqrt(betytwiss.ERRBETY[betytwiss.indx[name]]**2+betytwiss.STDBETY[betytwiss.indx[name]]**2)
			betmdl=betytwiss.BETYMDL[betytwiss.indx[name]]

			betp=normal_pro.BETY[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETY[errbetamax.indx[name]]-errbetamin.BETY[errbetamin.indx[name]])

			betb=back_pro.BETY[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETY[errbetamaxb.indx[name]]-errbetaminb.BETY[errbetaminb.indx[name]])
	
			databety.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
			
		# for ampbety
		print "filling table for ampbetay"
		for bpm in ampbpmsbety:
			name=bpm[1]
			s=ampbetytwiss.S[ampbetytwiss.indx[name]]
			#if ampbetytwiss.S[ampbetytwiss.indx[name]]<ampbetytwiss.S[ampbetytwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=ampbetytwiss.BETY[ampbetytwiss.indx[name]]
			errbet=ampbetytwiss.BETYSTD[ampbetytwiss.indx[name]]
			betmdl=betytwiss.BETYMDL[ampbetytwiss.indx[name]]

			betp=normal_pro.BETY[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETY[errbetamax.indx[name]]-errbetamin.BETY[errbetamin.indx[name]])

			betb=back_pro.BETY[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETY[errbetamaxb.indx[name]]-errbetaminb.BETY[errbetaminb.indx[name]])
	
			dataampbety.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		
		# for alfy
		print "filling table for alfy"
		for bpm in bpmsbety:
			name=bpm[1]
			s=betytwiss.S[betytwiss.indx[name]]
			#if betytwiss.S[betytwiss.indx[name]]<betytwiss.S[betytwiss.indx[element]]:s=s+twisstwiss.LENGTH
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
			#if Dx.S[Dx.indx[name]]<Dx.S[Dx.indx[element]]:s=s+twisstwiss.LENGTH
			dx=Dx.DX[Dx.indx[name]]
		        errdx=Dx.STDDX[Dx.indx[name]]
		        dxmdl=Dx.DXMDL[Dx.indx[name]]
			 
		        dpx=Dx.DPX[Dx.indx[name]]
		        dpxmdl=Dx.DPXMDL[Dx.indx[name]]
		        dpxerr="0"
			 
		        dxp=normal_pro.DX[normal_pro.indx[name]]
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
		        s=betxtwiss.S[Dy.indx[name]]
			#if Dy.S[Dy.indx[name]]<Dy.S[Dy.indx[element]]:s=s+twisstwiss.LENGTH
		        dx=Dy.DY[Dy.indx[name]]
		        errdx=Dy.STDDY[Dy.indx[name]]
		        dxmdl=Dy.DYMDL[Dy.indx[name]]
			 
		        dpx=Dy.DPY[Dy.indx[name]]
		        dpxmdl=Dy.DPYMDL[Dy.indx[name]]
		        dpxerr="0"
		        
		        dxp=normal_pro.DY[normal_pro.indx[name]]
		        errdxp=abs(errdmax.DY[errdmax.indx[name]]-errdmin.DY[errdmin.indx[name]])

		        dpxp=normal_pro.DPY[normal_pro.indx[name]]
		        dpxperr="0"
			 
		        dxb=back_pro.DY[back_pro.indx[name]]
		        errdxb=abs(errdmaxb.DY[errdmaxb.indx[name]]-errdminb.DY[errdminb.indx[name]])

		        dpxb=-back_pro.DPY[back_pro.indx[name]]
		        dpxberr="0"

			datady.append(name+" "+str(s)+" "+str(dx)+" "+str(errdx)+" "+str(dxmdl)+" "+str(dxp)+" "+str(errdxp)+" "+str(dxb)+" "+str(errdxb))
			datadpy.append(name+" "+str(s)+" "+str(dpx)+" "+str(dpxerr)+" "+str(dpxmdl)+" "+str(dpxp)+" "+str(dpxperr)+" "+str(dpxb)+" "+str(dpxberr))

		# for phasex total
		print "filing table for total horizontal phase"
		
		firstbpmx=bpmsphasex[0][1]
		endbpmx=bpmsphasex[len(bpmsphasex)-1][1]

	#	print len(bpmsphasex)
		#sys.exit()
	
		for bpms in bpmsphasex:
			bpm=bpms[1]
			S=filephasextot.S[filephasextot.indx[bpm]]
			prop=(normal_pro.MUX[normal_pro.indx[bpm]]-normal_pro.MUX[normal_pro.indx[firstbpmx]]) %1
			back=(1-(back_pro.MUX[back_pro.indx[bpm]]-back_pro.MUX[back_pro.indx[endbpmx]])) %1
			measured=(filephasextot.PHASEX[filephasextot.indx[bpm]]-filephasextot.PHASEX[filephasextot.indx[firstbpmx]])%1
			
			# difference experiment propogation
			MdiffP=(measured-prop ) %1
			if MdiffP >  0.5: MdiffP=  MdiffP -1
			elif MdiffP < -0.5: MdiffP=  MdiffP +1

			# difference experiment back propgation
			MdiffB=(measured-back) %1
			if MdiffB >  0.5: MdiffB=  MdiffB -1
			elif MdiffB < -0.5: MdiffB=  MdiffB +1
		
			errtotphasex=filephasextot.STDPHX[filephasextot.indx[bpm]]
			
			dataphasextot.append(bpm+' '+str(S)+' '+str(measured)+' '+str(errtotphasex)+' '+str(prop)+' '+str(0)+' '+str(MdiffP)+' '+str(back)+' '+str(0)+' '+str(MdiffB))
			# for phasey total
	
		firstbpmy=bpmsphasey[0][1]
		endbpmy=bpmsphasey[len(bpmsphasey)-1][1]
		
		print "filling table for total vertical phase"
		for bpms in bpmsphasey:
			bpm=bpms[1]
			S=filephaseytot.S[filephaseytot.indx[bpm]]
			prop=(normal_pro.MUY[normal_pro.indx[bpm]]-normal_pro.MUY[normal_pro.indx[firstbpmy]]) %1
			back=(1-(back_pro.MUY[back_pro.indx[bpm]]-back_pro.MUY[back_pro.indx[endbpmy]])) %1
			measured=(filephaseytot.PHASEY[filephaseytot.indx[bpm]]-filephaseytot.PHASEY[filephaseytot.indx[firstbpmy]])%1
			
			# difference experiment propogation
			MdiffP=(measured-prop ) %1
			if MdiffP >  0.5: MdiffP=  MdiffP -1
			elif MdiffP < -0.5: MdiffP=  MdiffP +1

			# difference experiment back propgation
			MdiffB=(measured-back) %1
			if MdiffB >  0.5: MdiffB=  MdiffB -1
			elif MdiffB < -0.5: MdiffB=  MdiffB +1

			errtotphasey=filephaseytot.STDPHY[filephaseytot.indx[bpm]]

			dataphaseytot.append(bpm+' '+str(S)+' '+str(measured)+' '+str(errtotphasey)+' '+str(prop)+' '+str(0)+' '+str(MdiffP)+' '+str(back)+' '+str(0)+' '+str(MdiffB))

		#for coupling
		print "filling table for coupling"
		for bpm in coupling[0]:
			#S=filecouple.S[filecouple.indx[bpm]]
			S=coupling[1][bpm][0]
			
			f1001=sqrt(coupling[1][bpm][1]**2+coupling[1][bpm][2]**2)
			f1001re=coupling[1][bpm][1]
			f1001im=coupling[1][bpm][2]
			stdf1001=abs(coupling[1][bpm][3])

			f1001p=abs(coupling[1][bpm][4])
			f1001pre=coupling[1][bpm][4].real
			f1001pim=coupling[1][bpm][4].imag
			stdf1001p=abs(coupling[1][bpm][5])

		        f1001b=abs(coupling[1][bpm][6])
		        f1001bre=coupling[1][bpm][6].real
			f1001bim=coupling[1][bpm][6].imag
		        stdf1001b=abs(coupling[1][bpm][7])

			f1001m=coupling[1][bpm][8]

			f1010=sqrt(coupling[1][bpm][9]**2+coupling[1][bpm][10]**2)
			f1010re=coupling[1][bpm][9]
			f1010im=coupling[1][bpm][10]
			stdf1010=abs(coupling[1][bpm][11])

			f1010m=coupling[1][bpm][12]
			
			f1010p=abs(coupling[1][bpm][13])
			f1010pre=coupling[1][bpm][13].real
			f1010pim=coupling[1][bpm][13].imag
			stdf1010p=coupling[1][bpm][14]
		
		
			f1010b=abs(coupling[1][bpm][15])
			f1010bre=coupling[1][bpm][15].real
			f1010bim=coupling[1][bpm][15].imag
			stdf1010b=abs(coupling[1][bpm][16])
		

			#print f1001pre,f1001bre,f1001re

			#print bpm,f1001re,f1001pre,f1001bre
			#sys.exit()
			
			coupledata.append(bpm+' '+str(S)+' '+str(f1001)+' '+str(stdf1001)+' '+str(f1001re)+' '+str(f1001im)+' '+str(f1001p)+' '+str(f1001pre)+' '+str(f1001pim)+' '+str(stdf1001p)+' '+str(f1001m)+' '+str(f1010)+' '+str(stdf1010)+' '+str(f1010re)+' '+str(f1010im)+' '+str(f1010p)+' '+str(f1010pre)+' '+str(f1010pim)+' '+str(stdf1010p)+' '+str(f1010m)+' '+str(f1001b)+' '+str(f1001bre)+' '+str(f1001bim)+' '+str(stdf1001b)+' '+str(f1010b)+' '+str(f1010bre)+' '+str(f1010bim)+' '+str(stdf1010b))
                 #########
			#sys.exit()
		flag=0


        #=> switch only for element - to - instrument
	elif passs==1:
        #else:

	        flag=1

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
	
		[hor,ver,dp]=getTwiss(savepath+"/twiss_"+str(namename)+".dat",namename)
		[horb,verb,dpb]=getTwiss(savepath+"/twiss_"+namename+"_back_rev.dat",namename)
		print "print the values"

		fileresul.write('@ POSSTART %le '+str(twisstwiss.S[twisstwiss.indx[startbpm]])+'\n')
		fileresul.write('@ BPMSTART %s '+str(startbpm)+'\n')
		fileresul.write('@ POSELEMENT %le '+str(twisstwiss.S[twisstwiss.indx[namename]])+'\n')
		fileresul.write('@ POSEND %le '+str(twisstwiss.S[twisstwiss.indx[startbpm]])+'\n')
		fileresul.write('@ BPMEND %s '+str(endbpm)+'\n')
		fileresul.write("@ INSTRUMENT  %s "+namename+"\n")
		fileresul.write("@ S %le  "+str(hor[2])+"\n")

		errbetxp=abs(errbetamax.BETX[errbetamax.indx[namename]]-errbetamin.BETX[errbetamin.indx[namename]])
		errbetxb=abs(errbetamaxb.BETX[errbetamaxb.indx[namename]]-errbetaminb.BETX[errbetaminb.indx[namename]])
		errbetyp=abs(errbetamax.BETY[errbetamax.indx[namename]]-errbetamin.BETY[errbetamin.indx[namename]])
		errbetyb=abs(errbetamaxb.BETY[errbetamaxb.indx[namename]]-errbetaminb.BETY[errbetaminb.indx[namename]])

	
		fileresul.write("* METHOD  BETX ERRBETX  ALFX   BETY  ERRBETY   ALFY   DX   DPX   DY   DPY F1001  F1001STD F1010  F1010STD\n")
		fileresul.write("$  %s  %le  %le %le  %le  %le  %le %le  %le  %le  %le %le  %le  %le  %le\n")

		if coupleswitch==0:
			fileresul.write("normal "+str(hor[0])+" "+str(errbetxp)+" "+str(hor[1])+" "+str(ver[0])+" "+str(errbetyp)+" "+str(ver[1])+" "+str(dp[0])+" "+str(dp[1])+" "+str(dp[2])+" "+str(dp[3])+" "+str(0)+" "+str(0)+" "+str(0)+" "+str(0)+"\n")
			fileresul.write("back "+str(horb[0])+" "+str(errbetxb)+" "+str(horb[1])+" "+str(verb[0])+" "+str(errbetyb)+" "+str(verb[1])+" "+str(dpb[0])+" "+str(dpb[1])+" "+str(dpb[2])+" "+str(dpb[3])+" "+str(0)+" "+str(0)+" "+str(0)+" "+str(0)+"\n")
		else:
			
			fileresul.write("normal "+str(hor[0])+" "+str(errbetxp)+" "+str(hor[1])+" "+str(ver[0])+" "+str(errbetyp)+" "+str(ver[1])+" "+str(dp[0])+" "+str(dp[1])+" "+str(dp[2])+" "+str(dp[3])+" "+str(0)+" "+str(0)+" "+str(0)+" "+str(0)+"\n")
				      #+str(coupling[1][name][5])+" "+str(coupling[1][name][7])+" "+str(coupling[1][name][6])+" "+str(coupling[1][name][8])+"\n")
			fileresul.write("back "+str(horb[0])+" "+str(errbetxb)+" "+str(horb[1])+" "+str(verb[0])+" "+str(errbetyb)+" "+str(verb[1])+" "+str(dpb[0])+" "+str(dpb[1])+" "+str(dpb[2])+" "+str(dpb[3])+" "+str(0)+" "+str(0)+" "+str(0)+" "+str(0)+"\n")
			
		fileresul.close()

			        #### creating tables
		#
		# => opening all the results
		#

		# opening data
		normal_pro=twiss(savepath+'/twiss_'+namename+'.dat')
		back_pro=twiss(savepath+'/twiss_'+namename+'_back_rev.dat')

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

		# for coupling
		if coupleswitch==1:
			S=0#coupling[1][name][0]
			f1001=0
			f1010=0
			stdf1001=0
			stdf1010=0

			f1001p=0#coupling[1][name][5]
			f1010p=0#coupling[1][name][6]
			stdf1001p=0#coupling[1][name][7]
			stdf1010p=0#coupling[1][name][8]

			f1001m=0#coupling[1][name][9]
			f1010m=0#coupling[1][name][10]			
			
			coupledata.append(name+' '+str(S)+' '+str(f1001)+' '+str(stdf1001)+' '+str(f1001p)+' '+str(stdf1001p)+' '+str(f1001m)+' '+str(f1010)+' '+str(stdf1010)+' '+str(f1010p)+' '+str(stdf1010p)+' '+str(f1010m))

                 #########
	else:
		flag=0
        print namename+" is Finished"


print "printing tables to file"
createTables("sbsbetx_"+namename+".out",savepath,columnnamesbetx,variablenamebetx,databetx,mainvariablebetx,mainvaluebetx)
createTables("sbsampbetx_"+namename+".out",savepath,columnnamesampbetx,variablenameampbetx,dataampbetx,mainvariableampbetx,mainvalueampbetx)
createTables("sbsalfx_"+namename+".out",savepath,columnnamesalfx,variablenamealfx,dataalfx,mainvariablealfx,mainvaluealfx)
createTables("sbsbety_"+namename+".out",savepath,columnnamesbety,variablenamebety,databety,mainvariablebety,mainvaluebety)
createTables("sbsampbety_"+namename+".out",savepath,columnnamesampbety,variablenameampbety,dataampbety,mainvariableampbety,mainvalueampbety)
createTables("sbsalfy_"+namename+".out",savepath,columnnamesalfy,variablenamealfy,dataalfy,mainvariablealfy,mainvaluealfy)
createTables("sbsdx_"+namename+".out",savepath,columnnames1dx,variablename1dx,datadx,mainvariabledx,mainvaluedx)
createTables("sbsdpx_"+namename+".out",savepath,columnnames2dx,variablename2dx,datadpx,mainvariabledx,mainvaluedx)
createTables("sbsdy_"+namename+".out",savepath,columnnames1dy,variablename1dy,datady,mainvariabledy,mainvaluedy)
createTables("sbsdpy_"+namename+".out",savepath,columnnames2dy,variablename2dy,datadpy,mainvariabledy,mainvaluedy)
createTables("sbsphasext_"+namename+".out",savepath,columnnamesphasextot,variablenamephasextot,dataphasextot,mainvariablephasextot,mainvaluephasextot)
createTables("sbsphaseyt_"+namename+".out",savepath,columnnamesphaseytot,variablenamephaseytot,dataphaseytot,mainvariablephaseytot,mainvaluephaseytot)
createTables("sbscouple_"+namename+".out",savepath,columnnamescouple,variablenamecouple,coupledata,mainvariablecouple,mainvaluecouple)








if flag==1:
	print "Making summary report for instruments"
	filee = os.listdir(savepath)
	files = filter(lambda x: 'resul_' in x    , filee)
	#files=files+filter(lambda x: 'resul_BW' in x    , filee)
	#files=files+filter(lambda x: 'resul_IP' in x    , filee)
	#files=files+filter(lambda x: 'resul_BPM' in x    , filee)

	#dateto=date.today()
	year=date.today().year
	month=date.today().month
	day=date.today().day
	
	resul=open(savepath+"/summary_instruments_"+options.accel+".tfs","w")
	resul.write('@ NAME %s "summary instruments"\n')
	resul.write('@ CREATED %s "'+str(year)+'-'+str(month)+'-'+str(day)+'"\n')
	resul.write('@ LABEL1 %s "OK->DATA_VALID"\n')
	resul.write('@ LABEL1 %s "CAUTION->SELEMENT-SBPM at some distance"\n')
	resul.write('@ LABEL1 %s "NOT_VALID->SELEMENT-SBPM to far"\n')
	
	if options.accel=="LHCB2":

		start="MKI.A5R8.B2"
		beam="B2"
		resul.write('@ ACCELERATOR %s "LHC"\n')
		resul.write('@ BEAM %s "'+beam+'"\n')
		resul.write('@ START-ELEMENT %s "'+start+'"\n')
        
        elif options.accel=="LHCB1":

		start="MSIA.EXIT.B1"   
		beam="B1"
		resul.write('@ ACCELERATOR %s "LHC"\n')
		resul.write('@ BEAM %s "'+beam+'"\n')
		resul.write('@ START-ELEMENT %s "'+start+'"\n')
	
	resul.write('* NAME  S	    BETX  BETERRX   BETXMDL  BETY  BETERRY  BETYMDL  LEFTBPM   RIGHTBPM   VALID\n')
	resul.write('$ %s    %le    %le   %le    %le      %le   %le   %le   %s  %s    %s\n')

	for filee in files:

		print savepath+"/"+filee
		one=twiss(savepath+"/"+filee)
		names=one.METHOD

		posstart=int(one.POSSTART)
		poselement=int(one.POSELEMENT)
		posend=int(one.POSEND)


		
		if(abs(poselement-posstart)<300 and abs(posend-poselement)<300):
			tag="OK"
		elif(abs(poselement-posstart)<1000 and abs(posend-poselement)<1000):
			tag="CAUTION"
		else:
			tag="NOT_VALID"
		

		betx=[];errx=[];bety=[];erry=[]

		try:
			betx.append(one.BETX[0])
			errx.append(1/one.ERRBETX[0])
			bety.append(one.BETY[0])
			erry.append(1/one.ERRBETY[0])
			betx.append(one.BETX[1])
			errx.append(1/one.ERRBETX[1])
			bety.append(one.BETY[1])
			erry.append(1/one.ERRBETY[1])


			betxx=(1/(errx[0]+errx[1]))*(errx[0]*betx[0]+errx[1]*betx[1])
			betyy=(1/(erry[0]+erry[1]))*(erry[0]*bety[0]+erry[1]*bety[1])
			errbetx=sqrt((1/errx[0])**2+(1/errx[1])**2)/2
			errbety=sqrt((1/erry[0])**2+(1/erry[1])**2)/2
		except: # will enter only here when error =0
			betxx=(one.BETX[0]+one.BETX[1])/2
			betyy=(one.BETY[0]+one.BETY[1])/2
			errbetx=0
			errbety=0
		name4colli=one.INSTRUMENT

		betxmdl=twisstwiss.BETX[twisstwiss.indx[name4colli]]
		betymdl=twisstwiss.BETY[twisstwiss.indx[name4colli]]
		S=twisstwiss.S[twisstwiss.indx[name4colli]]

		resul.write(name4colli+" "+str(S)+" "+str(betxx)+" "+str(errbetx)+" "+str(betxmdl)+" "+str(betyy)+" "+str(errbety)+" "+str(betymdl)+" "+one.BPMSTART+" "+one.BPMEND+" "+str(tag)+"\n")

	print "summary table created"
	resul.close()


print "sbs is finished"




		



	
