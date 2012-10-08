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
#        4 April 2010              -Introducing the modifiers.madx     
#                                   file, it should be in twisspath or an 
#                                   empty file will be created.             
#                                   The MAD mask file has also been 
#                                   updated
#  !=> SegementBySegment_0.29.py : - Update for IP2  beam1 and IP8 beam 2... problem with S coordinate
#  !=>  SegementBySegment_0.30.py : - New filtering... + extending it to RHIC (double plane bpms only, IPS)
#  !=>  SegementBySegment_0.31.py : - New and easier filtering



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
                help="basic twiss file, the modifiers.madx is assumed to be in the same direcory",
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
parser.add_option("-x", "--take", # take or create mad input, default should be 0 for creating
                help="take or create madx 0/1",
                metavar="mad", default="0", dest="madpass")
parser.add_option("-z", "--switch", # assumes that output is same as input
                help="switch one/two beam 1/2",
                metavar="switch", default="1", dest="switch")
parser.add_option("-c", "--cuts", # assumes that output is same as input
                help="cuts for beta,disp,coupling bbxm,bbym,bxe,bye,dxm,dym,dxe,dye,f1001e,f1010e",
                metavar="cuts", default="0.8,0.8,20,20,2,2,0.5,0.5,5,5", dest="cuts")
parser.add_option("-w", "--w", # Path to Chromaticity functions
                    help="Path to  chromaticity functions, by default this is skiped",
                    metavar="wpath", default="0", dest="wpath")

parser.add_option("-e", "--elementswitch", # assumes that output is same as input
                help="switching between segment(0) or element (1)",
                metavar="holyswitch", default="0", dest="holyswitch")


(options, args) = parser.parse_args()

####### some defs
def modelIntersect(expbpms, model):
	bpmsin=[]
        for bpm in expbpms:

		#print bpm
            
		try:
			check=model.indx[bpm.replace("-","_").upper()]
			#print "Found"
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


def filterbpm(betaxx,betayy,commonbpms):

	#
	# Easy function to find bpm and filter on beta>0,beterr>0 and err<beta
	#
	#
	goodbpm=[]
	for bp in range(len(commonbpms)):
		bpm=commonbpms[bp][1]
		betax=betaxx.BETX[betaxx.indx[bpm]]
		betxmdl=betaxx.BETXMDL[betaxx.indx[bpm]]
		betay=betayy.BETY[betayy.indx[bpm]]
		betymdl=betayy.BETYMDL[betayy.indx[bpm]]
		
		errbetax=betaxx.ERRBETX[betaxx.indx[bpm]]
		stdbetax=betaxx.STDBETX[betaxx.indx[bpm]]
		errbetay=betayy.ERRBETY[betayy.indx[bpm]]
		stdbetay=betayy.STDBETY[betayy.indx[bpm]]

		errx=sqrt(errbetax**2+stdbetax**2)
		erry=sqrt(errbetay**2+stdbetay**2)

	
		if (betax>0) and (betay>0) and (errx<betax) and (erry<betay) and (errx>0) and (erry>0):
			goodbpm.append(bpm)
			
		else:
			print "bpm ",bpm," didt not pass the cuts"
	
	
	
	return goodbpm


def getTwiss(filee,element):

    filedatax=twiss(filee)

    try:
	 filedatax.indx[element]
    except:
         print element," is not in file ",filee, " please check system exit"
	 sys.exit()

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

def run4mad(path,hor,ver,hore,vere,dp,dpe,startbpm,endbpm,name, fs, exppath,twissfile):


    ##
    #  Copy getfterms.py locally to be used by MADx
    cpath=options.bb
    if "RHIC" in options.accel:
	    os.system('cp '+cpath+'/SegmentBySegment/getfterms_RHIC.py'+' '+path+'/')
    else:
	    os.system('cp '+cpath+'/SegmentBySegment/getfterms_Dxy.py'+' '+path+'/')	    
    
    # Copy the modifiers.madx file locally to be used by MADx
    
    if os.path.isfile(twisspath+'/modifiers.madx'):
      os.system('cp '+twisspath+'/modifiers.madx'+' '+path+'/')
    else :   #If the modifiers file does not exist create empty file 
      os.system('touch '+path+'/modifiers.madx')
    
    
    if options.accel=="LHCB2":

        dire=-1
        start="MKI.A5R8.B2"
	#start="BPM.28L1.B2"
        beam="B2"
        
    elif options.accel=="LHCB1":

        dire=1
	start="MSIA.EXIT.B1"   #  compatible with repository
#	start="IP2"   #  compatible with repository
        beam="B1"

    elif options.accel=="RHICB":

        dire=1
	start="G6_BX"
	path4rhic=options.bb+"/MODEL/"+options.accel
	path4asc=os.path.dirname(options.twiss)
	beam="blue.asc"


    elif options.accel=="RHICY":

        dire=-1
	start="G6_BX"
	path4rhic=options.bb+"/MODEL/"+options.accel
	path4asc=os.path.dirname(options.twiss)
	beam="yellow.asc"
   

    ### Check chromatic functions in wpath
    if options.wpath!="0":
      print "Chrom path,", options.wpath
      try:
        wx=twiss(options.wpath+"/wx.out")   
        wy=twiss(options.wpath+"/wy.out")
      except:
        print "No Chromatic functions (wx,wy) available at ",options.wpath
        options.wpath="0" # This means set wx,wy=0,0
        
    if options.wpath!="0":
      try:
        wxs=wx.WX[wx.indx[startbpm]]
        phixs=wx.PHIX[wx.indx[startbpm]]
      except:
        print "Start BPM, ",startbpm," not in WX file"
        wxs=0
        phixs=0

      try:
        wys=wy.WY[wy.indx[startbpm]]
        phiys=wy.PHIY[wy.indx[startbpm]]
      except:
        print "Start BPM, ",startbpm," not in WY file"
        wys=0
        phiys=0
    else:
       wxs=0
       phixs=0
       wys=0
       phiys=0

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
    file4nad.write('    -e \'s/%STARTFROM/\''+str(startbpm.replace("-","_"))+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ENDAT/\''+str(endbpm.replace("-","_"))+'\'/g\' \\\n')
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
    file4nad.write('    -e \'s/%WX/\''+str(wxs)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PHIX/\''+str(phixs)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%WY/\''+str(wys)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PHIY/\''+str(phiys)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%WPATH/\'\"'+str(options.wpath.replace('/','\/'))+'\"\'/g\' \\\n')
    


    if "RHIC" in options.accel:
	    #print twissfile
            file4nad.write('    -e \'s/%QX/\''+str(twissfile.Q1)+'\'/g\' \\\n')
            file4nad.write('    -e \'s/%QY/\''+str(twissfile.Q2)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%CHROMX/\''+str(twissfile.DQ1)+'\'/g\' \\\n')
            file4nad.write('    -e \'s/%CHROMY/\''+str(twissfile.DQ2)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%BBPATH/\'\"'+str(path4rhic.replace('/','\/'))+'\"\'/g\' \\\n')
	    file4nad.write('    -e \'s/%MPATH/\'\"'+str(path4asc.replace('/','\/'))+'\"\'/g\' \\\n')
	    file4nad.write('    -e \'s/%MASS/\''+str(twissfile.MASS)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%ENERGY/\''+str(twissfile.ENERGY)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%GAMMA/\''+str(twissfile.GAMMA)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%CHARGE/\''+str(twissfile.CHARGE)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%EMX/\''+str(twissfile.EX)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%EMY/\''+str(twissfile.EY)+'\'/g\' \\\n')
	    file4nad.write('    -e \'s/%SIGE/\''+str(twissfile.SIGE)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%EXP/\'\"'+exppath.replace("/","\/")+'\"\'/g\' \\\n')
    
    
    
    if "RHIC" in options.accel:

	   file4nad.write('<'+cpath+'/SegmentBySegment/'+'/job.InterpolateBetas.0_1_RHIC.mask > '+path+'/t_'+str(name)+'.madx \n') 

    else:
    
	    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/job.InterpolateBetas.0_2.mask > '+path+'/t_'+str(name)+'.madx \n')

    file4nad.close()
    
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))
    runmad(path,name)
    #sys.exit()

   
    
def runmad(path,name):
	
	os.system(options.mad+'madx < '+path+'t_'+str(name)+'.madx')
	
   
  
def run4plot(path,spos,epos,beta4plot,cpath,meapath,name,qx,qy,accel):

    filename=path+'/var4plot.sh'
    file4nad=open(filename,'w')
    file4nad.write('sed -e \'s/%PATH/\'\"'+str(path.replace("/","\/"))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%EndPoint/\''+str(epos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%StartPoint/\''+str(spos)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LABEL/\''+str(name)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%ACCEL/\''+str(options.accel)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%BETA/\''+str(beta4plot)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%QX/\''+str(qx)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%QY/\''+str(qy)+'\'/g\' \\\n')    
    file4nad.write('    -e \'s/%MEA/\'\"'+str(meapath.replace("/","\/"))+'\"\'/g\' \\\n')
    if (name=="IP8" and accel=="LHCB2") or (name=="IP2" and accel=="LHCB1"):
	    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.IP2IP8.0_2.mask > '+path+'/gplot_'+name)
    elif "RHIC" in options.accel:
	     file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.0_1_RHIC.mask > '+path+'/gplot_'+name)
    else:
	    file4nad.write('<'+cpath+'/SegmentBySegment/'+'/gplot.0_2.mask > '+path+'/gplot_'+name)
	    

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


path=options.path
beamswitch=options.switch
##########
#
# LOADING DATA
#
##########

filedatax=twiss(path+"/getbetax.out")
betxtwiss=filedatax
QXX=filedatax.Q1
QYY=filedatax.Q2
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
		
### check if dispersion exist
try:
	filedx=twiss(path+"/getDx.out")
	filedy=twiss(path+"/getDy.out")
	disp=1
        print "Disp files OK", filedx.DX[0], filedx.NAME[0]
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

# checking if end element is BPM or instrument (IP,collimators)
elementswitch=-1
savepath=options.SAVE
list2run=options.segf.split(',')
cuts=options.cuts.split(',')

start={}
end={}
names=[]


twissfile=options.twiss
twisspath=os.path.dirname(twissfile)+'/'
if twissfile=="./":
	twissfile=options.bb+"/MODEL/"+options.accel+"/nominal.opt/twiss.dat"

twisstwiss=twiss(twissfile)



	
#if options.holyswitch=="0":
if (len(list2run)/3)!=1:
	print "You didnt gave the right input for segment start,end,name => System exit"
	sys.exit()
for run in range(len(list2run)-2):
		
	start[list2run[run+2]]=list2run[run]
	end[list2run[run+2]]=list2run[run+1]
	names.append(list2run[run+2])


#else:
	#for run in range(len(list2run)):
		#start[list2run[run]]=twisstwiss.NAME[0]
		#end[list2run[run]]=twisstwiss.NAME[len(twisstwiss.NAME)-1]
		#names.append(list2run[run])


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



######## big loop

####### Difference between segment or instrument #######
# For segment => Startbpm,endbpm,NameOfSegment
# For instrument => StartBPM,endBpm, instrument (must be excisting name)
####### 

for namename in names:



	passs=1
	commonbpms=intersect([betxtwiss,betytwiss,twisstwiss])
	bpms=[]
	
	for i in range(len(commonbpms)):
		bpms.append(commonbpms[i][1])

	#print commonbpms
	#print start[namename],end[namename]


 	
	if options.holyswitch=="0":
		elementswitch=0
		print "Segment has been choosen"

		
		
		#print start[namename] not in bpms,end[namename] not in bpms
		if (start[namename] not in bpms) or (end[namename] not in bpms) :
				print "The BPMS choosen are not in the measurement files ... please change borders"
				if start[namename] not in bpms:
					print "startbpm is not in measurement",start[namename]
				if end[namename] not in bpms:
					print "endbpm is not in measurement",	end[namename]				
				sys.exit()


		#finding good bpm
		goodbpms=filterbpm(betxtwiss,betytwiss,commonbpms)
		
		if (start[namename] not in goodbpms) or (end[namename] not in goodbpms):
			print "Not enough bpms for propogation...  please increase border => BPMS lost due to cut"
			if start[namename] not in bpms:
					print "startbpm did not pass the cuts",start[namename]
			if end[namename] not in bpms:
					print "endbpm dit not pass the cuts",end[namename]
					
			print "Skipping ",namename
			passs=0
			sys.exit()
		else:
			startbpm=start[namename]
			endbpm=end[namename]
			#print startbpm,endbpm
			#print goodbpms
			

	else:
		elementswitch=1
		goodbpms=filterbpm(betxtwiss,betytwiss,commonbpms)
		if (start[namename] not in goodbpms) or (end[namename] not in goodbpms):
			print "Not enough bpms for propogation...  please increase border => BPMS lost due to cut"
			if start[namename] not in bpms:
					print "startbpm did not pass the cuts",start[namename]
			if end[namename] not in bpms:
					print "endbpm dit not pass the cuts",end[namename]
					
			print "Skipping ",namename
			passs=0
			sys.exit()
		else:
			startbpm=start[namename]
			endbpm=end[namename]
			#print startbpm,endbpm
			#print goodbpms



	


	#gathering data
	hor=[betxtwiss.BETX[betxtwiss.indx[startbpm]],
	     betxtwiss.STDBETX[betxtwiss.indx[startbpm]],
	     betxtwiss.ALFX[betxtwiss.indx[startbpm]],
	     betxtwiss.ERRALFX[betxtwiss.indx[startbpm]]]
	ver=[betytwiss.BETY[betytwiss.indx[startbpm]],
	     betytwiss.STDBETY[betytwiss.indx[startbpm]],
	     betytwiss.ALFY[betytwiss.indx[startbpm]],
	     betytwiss.ERRALFY[betytwiss.indx[startbpm]]]
	hore=[betxtwiss.BETX[betxtwiss.indx[endbpm]],
	     betxtwiss.STDBETX[betxtwiss.indx[endbpm]],
	     betxtwiss.ALFX[betxtwiss.indx[endbpm]],
	     betxtwiss.ERRALFX[betxtwiss.indx[endbpm]]]
	vere=[betytwiss.BETY[betytwiss.indx[endbpm]],
	     betytwiss.STDBETY[betytwiss.indx[endbpm]],
	     betytwiss.ALFY[betytwiss.indx[endbpm]],
	     betytwiss.ERRALFY[betytwiss.indx[endbpm]]]



			

	if disp==1:

		try: #horizontal
			filedx.indx[startbpm]
			dxpass=1
		except:
			print "startbpm for horizontal dispersion not found"
			dxpass=0

		try: #vertical
			filedy.indx[startbpm]
			dypass=1
		except:
			print "startbpm for  vertical dispersion not found"
			dypass=0

		try: #horizontal
			filedx.indx[endbpm]
			dxpasse=1
		except:
			print "endbpm for horizontal dispersion not found"
			dxpasse=0

		try: #vertical
			filedy.indx[endbpm]
			dypasse=1
		except:
			print "endbpm for  vertical dispersion not found"
			dypasse=0

		if dxpass==1:
			dxx=filedx.DX[filedx.indx[startbpm]]
			dxp=filedx.DPX[filedx.indx[startbpm]]
			dxe=filedx.STDDX[filedx.indx[startbpm]]
		else:
			dxx=0
			dxp=0
			dxe=0

		if dypass==1:
			dyy=filedy.DY[filedy.indx[startbpm]]
			dyp=filedy.DPY[filedy.indx[startbpm]]
			dye=filedy.STDDY[filedy.indx[startbpm]]
		else:
			dyy=0
			dyp=0
			dye=0			
		
		dp=[dxx,
		    dxp,
		    dyy,
		    dyp,
		    dxe,
		    dye]

		if dxpass==1:
			dxx=filedx.DX[filedx.indx[endbpm]]
			dxp=filedx.DPX[filedx.indx[endbpm]]
			dxe=filedx.STDDX[filedx.indx[endbpm]]
		else:
			dxx=0
			dxp=0
			dxe=0

		if dypass==1:
			dyy=filedy.DY[filedy.indx[endbpm]]
			dyp=filedy.DPY[filedy.indx[endbpm]]
			dye=filedy.STDDY[filedy.indx[endbpm]]
		else:
			dyy=0
			dyp=0
			dye=0	
		    
		dpe=[dxx,
		    dxp,
		    dyy,
		    dyp,
		    dxe,
		    dye]		   
		
	else:
		dp=[0,0,0,0,0,0]
		dpe=[0,0,0,0,0,0]
		print "No dispersion"




		
	
	


	if passs==1:

		if coupleswitch==1:
			print "coupling,"
			coupling=[[],[]]
			try:
				filecouple.indx[startbpm]
				cpass=1
			except:
				cpass=0
				print "startbpm ",startbpm," not found in coupling measurement => values=0"
				
				
			if cpass==1:
				f1001r=filecouple.F1001R[filecouple.indx[startbpm]]
				f1001i=filecouple.F1001I[filecouple.indx[startbpm]]
				f1010r=filecouple.F1010R[filecouple.indx[startbpm]]
				f1010i=filecouple.F1010I[filecouple.indx[startbpm]]
			else:
				f1001r=0
				f1001i=0
				f1010r=0
				f1010i=0
			fs=[f1001r,f1001i,f1010r,f1010i]

	        else:
			coupling=[[],[]]
			fs=[0,0,0,0]
			print "No coupling"

		#sys.exit()
                print "madpass", options.madpass		
		if str(options.madpass)=="0":
                        print "Going to run4mad"
			run4mad(savepath,hor,ver,hore,vere,dp,dpe,startbpm,endbpm,namename, fs, options.path+"/",twisstwiss)

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
		#print savepath+"/twiss_"+namename+".dat"
		m=twiss(savepath+"/twiss_"+namename+".dat")
		#sys.exit()

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
		[hor,ver,dp]=getTwiss(twissfile,startbpm)
		beta4plot=hor[2]

	        ######################
		if options.gra=='1':
			
			run4plot(savepath,startpos,endpos,beta4plot,options.bb,path,namename,QXX,QYY,options.accel)


	

	        #### creating tables
		#
		# => opening all the results
		#


		# results from propagation
		normal_pro=twiss(savepath+'/twiss_'+namename+'.dat')
		back_pro=twiss(savepath+'/twiss_'+namename+'_back_rev.dat')
		bpmsbetx=modelIntersect(filedatax.NAME,back_pro)
		#print len(filedatax.NAME),len(back_pro.NAME),len(bpmsbetx)
		ampbpmsbetx=modelIntersect(ampbetxtwiss.NAME,back_pro)
		bpmsbety=modelIntersect(filedatay.NAME,back_pro)
		ampbpmsbety=modelIntersect(ampbetytwiss.NAME,back_pro)
		bpmsphasex=modelIntersect(filephasextot.NAME,back_pro)
		bpmsphasey=modelIntersect(filephaseytot.NAME,back_pro)
	
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
		#print  bpmsbetx
		#print normal_pro.NAME
		for bpm in bpmsbetx:
			name=bpm
			s=betxtwiss.S[betxtwiss.indx[name]]  

			bet=betxtwiss.BETX[betxtwiss.indx[name]]
			errbet=sqrt(betxtwiss.ERRBETX[betxtwiss.indx[name]]**2+betxtwiss.STDBETX[betxtwiss.indx[name]]**2)
			betmdl=betxtwiss.BETXMDL[betxtwiss.indx[name]]

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function
				name=name.replace("-","_")
				

			betp=normal_pro.BETX[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETX[errbetamax.indx[name]]-errbetamin.BETX[errbetamin.indx[name]])

			betb=back_pro.BETX[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETX[errbetamaxb.indx[name]]-errbetaminb.BETX[errbetaminb.indx[name]])
	
			databetx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		# for ampbetx
		
		print "filling table for ampbetx"
		for bpm in ampbpmsbetx:
			name=bpm
			s=ampbetxtwiss.S[ampbetxtwiss.indx[name]]
			#if ampbetxtwiss.S[ampbetxtwiss.indx[name]]<ampbetxtwiss.S[ampbetxtwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=ampbetxtwiss.BETX[ampbetxtwiss.indx[name]]
			errbet=ampbetxtwiss.BETXSTD[ampbetxtwiss.indx[name]]
			betmdl=ampbetxtwiss.BETXMDL[ampbetxtwiss.indx[name]]

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function
				name=name.replace("-","_")

			betp=normal_pro.BETX[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETX[errbetamax.indx[name]]-errbetamin.BETX[errbetamin.indx[name]])

			betb=back_pro.BETX[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETX[errbetamaxb.indx[name]]-errbetaminb.BETX[errbetaminb.indx[name]])
	
			dataampbetx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		
		# for alphax
		print "filling table for alfx"
		for bpm in bpmsbetx:
			name=bpm
			s=betxtwiss.S[betxtwiss.indx[name]]
			#if betxtwiss.S[betxtwiss.indx[name]]<betxtwiss.S[betxtwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=betxtwiss.ALFX[betxtwiss.indx[name]]
			errbet=sqrt(betxtwiss.ERRALFX[betxtwiss.indx[name]]**2+betxtwiss.STDALFX[betxtwiss.indx[name]]**2)
			betmdl=betxtwiss.ALFXMDL[betxtwiss.indx[name]]

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function
				name=name.replace("-","_")
			
			betp=normal_pro.ALFX[normal_pro.indx[name]]
			errbetp=abs(erralfmax.ALFX[erralfmax.indx[name]]-erralfmin.ALFX[erralfmin.indx[name]])

			betb=-back_pro.ALFX[back_pro.indx[name]]
			errbetb=abs(erralfmaxb.ALFX[erralfmaxb.indx[name]]-erralfminb.ALFX[erralfminb.indx[name]])
	
			dataalfx.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
					
                #for bety
		print "filling table for bety"
		for bpm in bpmsbety:
			name=bpm
			s=betytwiss.S[betytwiss.indx[name]]
			#if betytwiss.S[betytwiss.indx[name]]<betytwiss.S[betytwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=betytwiss.BETY[betytwiss.indx[name]]
			errbet=sqrt(betytwiss.ERRBETY[betytwiss.indx[name]]**2+betytwiss.STDBETY[betytwiss.indx[name]]**2)
			betmdl=betytwiss.BETYMDL[betytwiss.indx[name]]

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function
				name=name.replace("-","_")

			betp=normal_pro.BETY[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETY[errbetamax.indx[name]]-errbetamin.BETY[errbetamin.indx[name]])

			betb=back_pro.BETY[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETY[errbetamaxb.indx[name]]-errbetaminb.BETY[errbetaminb.indx[name]])
	
			databety.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
			
		# for ampbety
		print "filling table for ampbetay"
		for bpm in ampbpmsbety:
			name=bpm
			s=ampbetytwiss.S[ampbetytwiss.indx[name]]
			#if ampbetytwiss.S[ampbetytwiss.indx[name]]<ampbetytwiss.S[ampbetytwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=ampbetytwiss.BETY[ampbetytwiss.indx[name]]
			errbet=ampbetytwiss.BETYSTD[ampbetytwiss.indx[name]]
			betmdl=betytwiss.BETYMDL[ampbetytwiss.indx[name]]

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function like i want:-)
				name=name.replace("-","_")

			betp=normal_pro.BETY[normal_pro.indx[name]]
			errbetp=abs(errbetamax.BETY[errbetamax.indx[name]]-errbetamin.BETY[errbetamin.indx[name]])

			betb=back_pro.BETY[back_pro.indx[name]]
			errbetb=abs(errbetamaxb.BETY[errbetamaxb.indx[name]]-errbetaminb.BETY[errbetaminb.indx[name]])
	
			dataampbety.append(name+" "+str(s)+" "+str(bet)+" "+str(errbet)+" "+str(betmdl)+" "+str(betp)+" "+str(errbetp)+" "+str(betb)+" "+str(errbetb))
		
		# for alfy
		print "filling table for alfy"
		for bpm in bpmsbety:
			name=bpm
			s=betytwiss.S[betytwiss.indx[name]]
			#if betytwiss.S[betytwiss.indx[name]]<betytwiss.S[betytwiss.indx[element]]:s=s+twisstwiss.LENGTH
			bet=betytwiss.ALFY[betytwiss.indx[name]]
			errbet=sqrt(betytwiss.ERRALFY[betytwiss.indx[name]]**2+betytwiss.STDALFY[betytwiss.indx[name]]**2)
			betmdl=betytwiss.ALFYMDL[betytwiss.indx[name]]

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function
				name=name.replace("-","_")

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

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function
				name=name.replace("-","_")
			 
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
		        s=betytwiss.S[Dy.indx[name]]
			#if Dy.S[Dy.indx[name]]<Dy.S[Dy.indx[element]]:s=s+twisstwiss.LENGTH
		        dx=Dy.DY[Dy.indx[name]]
		        errdx=Dy.STDDY[Dy.indx[name]]
		        dxmdl=Dy.DYMDL[Dy.indx[name]]
			 
		        dpx=Dy.DPY[Dy.indx[name]]
		        dpxmdl=Dy.DPYMDL[Dy.indx[name]]
		        dpxerr="0"

			if "RHIC" in options.accel: # not nice i know :-(, dict doesn't function
				name=name.replace("-","_")
		        
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
		coupledata=[]
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








if flag==1:
	print "Making summary report for instruments"
	filee = os.listdir(savepath)
	files = filter(lambda x: 'resul_' in x    , filee)

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
		resul.write('@ STARTELEMENT %s "'+start+'"\n')
        
        elif options.accel=="LHCB1":

		start="MSIA.EXIT.B1"   
		beam="B1"
		resul.write('@ ACCELERATOR %s "LHC"\n')
		resul.write('@ BEAM %s "'+beam+'"\n')
		resul.write('@ STARTELEMENT %s "'+start+'"\n')

	elif  options.accel=="RHICB":

		start="G6_BX"   
		beam="BLUE"
		resul.write('@ ACCELERATOR %s "LHC"\n')
		resul.write('@ BEAM %s "'+beam+'"\n')
		resul.write('@ STARTELEMENT %s "'+start+'"\n')

	elif  options.accel=="RHICY":

		start="G6_BX"   
		beam="YELLOW"
		resul.write('@ ACCELERATOR %s "LHC"\n')
		resul.write('@ BEAM %s "'+beam+'"\n')
		resul.write('@ STARTELEMENT %s "'+start+'"\n')	
	
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




		



	
