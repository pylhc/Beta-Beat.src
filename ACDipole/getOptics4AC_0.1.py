################################################################
#                                                              #
#  @ Glenn Vanbavinckhove  (gvanbavi@cern.ch)=> Date: 11/02/10 #
#  @ Ryoichi Miyamoto (miyamoto@bnl.gov)                       #
#                                                              #
################################################################
#
#  !=> getOptics4AC_0.0.py : - Construction of the program (11/02/10)
#  !=> getOptics4AC_0.1.py : - Changing starting bpm for ac dipole (23/02/10)

# imports
from metaclass import twiss
from optparse import OptionParser
from math import *
import sys,os

# option parser
parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                metavar="ACCEL", default="LHCB1",dest="accel")
parser.add_option("-p", "--path",
                help="Path to output files of GetLLM",
                metavar="PATH", default="./",dest="path")
parser.add_option("-b", "--bbsrouce",
                help="beta beat source",
                metavar="bb", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/", dest="bb")
parser.add_option("-o", "--output path",
                help="Where to put the output files",
                metavar="output", default="./", dest="output")
parser.add_option("-d", "--drive tune",
                help="input of drive tune for horizontal and vertical",
                metavar="dtune", default="0.28,0.31", dest="dtune")
parser.add_option("-m", "--madx",
		help="madx executable",
		metavar="madx_var",default="/afs/cern.ch/group/si/slap/bin/madx_dev",dest="madx_var")

(options, args) = parser.parse_args()

#function
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


def getbeta(lambdad,phasebpm,phaseac,phasedriven,driventune,drivenbeta,sign):

    #print "phase ac ",phaseac%(2*pi)
    if sign=="+":
        #print phasebpm-phaseac+phasedriven+pi*driventune
        forcos=2*(phasebpm-phaseac+phasedriven+pi*driventune)
    else:
        forcos=2*(phasebpm-phaseac+phasedriven-pi*driventune)

    up=1+lambdad**2+2*lambdad*cos(forcos)
    down=1-lambdad**2

    
    
    beta=(up/down)*drivenbeta

    #print drivenbeta,up,down,beta,2*lambdad*cos(forcos),lambdad
    
    return beta

def getalfa(lambdad,phasebpm,phaseac,phasedriven,driventune,drivenalfa,sign):

    if sign=="+":
        forcos=2*(phasebpm-phaseac+phasedriven+pi*driventune)
    else:
        forcos=2*(phasebpm-phaseac+phasedriven-pi*driventune)

    up1=1+lambdad**2+2*lambdad*cos(forcos)
    up2=2*lambdad*sin(forcos)
    down=1-lambdad**2

    alfa=(up1/down)*drivenalfa+(up2/down)

    return alfa


def getphase(sign,lambdad,phasebpm,phaseac,phasedriven,driventune):

    if sign=="+":
        if phasebpm<pi:
            up=(1-lambdad**2)*sin(phasebpm)
            down=(1+lambdad**2)*cos(phasebpm)+2*lambdad*cos(phasebpm-2*phaseac+2*phasedriven+2*pi*driventune)
            phase=atan(up/down)
        else:
            up=(1-lambdad**2)*sin(phasebpm)
            down=(1+lambdad**2)*cos(phasebpm)+2*lambdad*cos(phasebpm-2*phaseac+2*phasedriven+2*pi*driventune)
            phase=atan(up/down)+pi*(phasebpm-pi)
    else:
        if phasebpm<pi:
            up=(1-lambdad**2)*sin(phasebpm-2*pi*driventune)
            down=(1+lambdad**2)*cos(phasebpm-2*pi*driventune)+2*lambdad*cos(phasebpm-2*phaseac+2*phasedriven)
            phase=atan(up/down)
        else:
            up=(1-lambdad**2)*sin(phasebpm-2*pi*driventune)
            down=(1+lambdad**2)*cos(phasebpm-2*pi*driventune)+2*lambdad*cos(phasebpm-2*phaseac+2*phasedriven)
            phase=atan(up/down)+pi*(phasebpm-2*pi*driventune-pi)+2*pi*driventune

    return phase


def getdrivenphase(phasemodel,tune,lambdad,driventune,Q):

    #phasetune=(phasemodel-pi*tune)
    phasetune=(phasemodel-pi*tune)%(2*pi)
    #print "phasetune",phasetune
    #print phasetune,"phasetune"
    lambdafun=(1+lambdad)/(1-lambdad)
    #print phasemodel*2*pi,lambdafun,lambdad
    #print phasetune>pi, phasetune<2*pi
    if phasetune<pi:
        #print "smaller ",atan(lambdafun*tan(phasetune))%pi
        phasedriven=(((atan(lambdafun*tan(phasetune)))%pi)+pi*driventune)%(2*pi)
    else:
        #print "other"
        #print "smaller ",atan(lambdafun*tan(phasetune))%pi
        phasedriven=(((atan(lambdafun*tan(phasetune)))%pi)+pi+pi*driventune)%(2*pi)

    print "driven phase",phasedriven

    return phasedriven
    


def getlambda(drivetune,tune):

    

    up=sin(pi*(drivetune-tune))
    down=sin(pi*(drivetune+tune))
    lambdad=up/down

    return lambdad

def getphasesum(data):

	names=data.NAME
	phasesum={}
	phasesum[names[0]]=0
	phasex=0
	phasey=0
	for bpm in range(len(names)-1):

		bpmname=names[bpm+1]
		try:
			phasex=phasex+data.PHASEX[data.indx[bpmname]]
#			print data.PHASEX[data.indx[bpmname]]
                        print phasex%(2*pi),bpmname,phasex
			
			phasesum[bpmname]=phasex
			
		except:
			phasey=phasey+data.PHASEY[data.indx[bpmname]]
			phasesum[bpmname]=phasey
	#sys.exit()	
			
	return phasesum	


def callandrunmadx(beam,plane,path,bsource,freetune,driventune,madx):

	maindir=os.path.dirname(path)
	acpath=maindir+"/ACDIPOLE.core/"

	if not os.path.exists(acpath):
		os.mkdir(acpath)

	if plane=="H":
		tunes=open(acpath+'/tunes'+beam+'_H.sh','w')
		shtunes=acpath+'/tunes'+beam+'_H.sh'
		madxpath=bsource+'/MODEL/'+beam+'/AC.Dipole/twiss.'+beam+'.hac.mask'
		madxpro='twiss.'+beam+'.hac.madx'
		print madxpath
		#if not os.path.exists(acpath+'/twiss.'+beam+'.hac.madx'):
			#os.system('cp '+madxpath+' '+acpath)
	else:
		tunes=open(acpath+'/tunes'+beam+'_V.sh','w')
		shtunes=acpath+'/tunes'+beam+'_V.sh'
		madxpath=bsource+'/MODEL/'+beam+'/AC.Dipole/twiss.'+beam+'.vac.mask'
		madxpro='twiss.'+beam+'.vac.madx'
		#if not os.path.exists(acpath+'/twiss.'+beam+'.vac.madx'):
			#os.system('cp '+madxpath+' '+acpath)




	print >> tunes,'sed -e \'s/%Q/\''+str(freetune)+'\'/g\' \\'
	print >> tunes,'    -e \'s/%D/\''+str(driventune)+'\'/g\' \\'
	print >> tunes,'    -e \'s/%PATH/\''+'\"'+str(acpath.replace('/','\/'))+'\"'+'\'/g\' \\'
	if plane=="H":
		print >> tunes,'<',madxpath,'> ',acpath+'/twiss.'+beam+'.hac.madx'
	else:
		print >> tunes,'<',madxpath,'> ',acpath+'/twiss.'+beam+'.vac.madx'		
	tunes.close()
	print 'chmod 777 '+shtunes
	os.system('chmod 777 '+shtunes)
	os.system(shtunes)

	print "Triggering madx job"
	print madx+'<'+acpath
	os.system(madx+'<'+acpath+'/'+madxpro)

	if plane=="H":
		twissfile=acpath+"/twiss."+beam+".hac.dat"
	else:
		twissfile=acpath+"/twiss."+beam+".vac.dat"

	twissac=twiss(twissfile)

	return twissac


def dosub(plane,real,model):
	
	bpms=intersect([real,model])
	data={}
	for bpm in bpms:

		s=bpm[0]
		bpm=bpm[1]

		if plane=="H":
			beta_real=real.BETX[real.indx[bpm]]
			
			alfa_real=real.ALFX[real.indx[bpm]]
			beta_model=model.BETX[model.indx[bpm]]
			alfa_model=model.ALFX[model.indx[bpm]]
			beta=(beta_real-beta_model)/beta_model
			#print beta
			alfa=alfa_real-alfa_model
			phase=0
		
		else:
			beta_real=real.BETY[real.indx[bpm]]
			alfa_real=real.ALFY[real.indx[bpm]]
			beta_model=model.BETY[model.indx[bpm]]
			alfa_model=model.ALFY[model.indx[bpm]]
			beta=(beta_real-beta_model)/beta_model
			alfa=alfa_real-alfa_model
			phase=0

		data[bpm]=[bpm,s,beta,alfa,phase]
	
	return [bpms,data]

	


# reading files
phasex=twiss(options.path+'/getphasex.out')
drivetunex=phasex.Q1
drivetuney=phasex.Q2
phasey=twiss(options.path+'/getphasey.out')
#phasext=twiss(options.path+'/getphaserx.out')
#phaseyt=twiss(options.path+'/getphasery.out')
betax=twiss(options.path+'/getbetax.out')
betay=twiss(options.path+'/getbetay.out')
ampbetax=twiss(options.path+'/getampbetax.out')
ampbetay=twiss(options.path+'/getampbetay.out')

#model=twiss(options.bb+'/MODEL/'+options.accel+'/nominal.opt/twiss_ac.dat')
allfiles=[phasex,phasey,betax,betay,ampbetax,ampbetay]

# writing files
betaxac=open(options.output+'/getbetax_acfree.out','w')
betayac=open(options.output+'/getbetay_acfree.out','w')

betaxac_mo=open(options.output+'/getbetax_acfree_mo.out','w')
betayac_mo=open(options.output+'/getbetay_acfree_mo.out','w')

print >> betaxac, "* NAME S BETAX ALFAX PHASEX BETAXM ALFAXM PHASEXM PHASED"
print >> betaxac, "$ %s   %le %le %le %le %le %le %le %le"

print >> betayac, "* NAME S BETAY ALFAY PHASEY BETAYM ALFAYM PHASEYM PHASED"
print >> betayac, "$ %s   %le %le %le %le %le %le %le %le"

print >> betaxac_mo, "* NAME S BETAX ALFAX PHASEX BETAXM ALFAXM PHASEXM PHASED"
print >> betaxac_mo, "$ %s   %le %le %le %le %le %le %le %le"

print >> betayac_mo, "* NAME S BETAY ALFAY PHASEY BETAYM ALFAYM PHASEYM PHASED"
print >> betayac_mo, "$ %s   %le %le %le %le %le %le %le %le"




# main part
bpms=intersect(allfiles)
#bpms=modelIntersect(bpms, model)
#nbpm="BPMYA.5L4.B1"
#phasextot=getphasesum(phasex)
#getphase_x=open('getphase4ryo.out','w')
#print len(bpms)
#print >>getphase_x , "NAME S PHASEXT PHASEXT_2PI"
#print >>getphase_x , "%s   %le %le   %le "
#sys.exit()
#for i in range(len(bpms)):
#	bpm=bpms[i][1]
#	print >>getphase_x ,bpms[i][1],bpms[i][0],phasextot[bpm],phasextot[bpm]%(2*pi)
#	print bpm,phasextot[bpm],bpms[0][1]

#getphase_x.close()
	
#phaseytot=getphasesum(phasey)

Qx=float(options.dtune.split(',')[0])
Qy=float(options.dtune.split(',')[1])
#lambdadx=getlambda(drivetunex,Qx)
#lambdady=getlambda(drivetuney,Qy)
#phasemodelx=model.MUX[model.indx[nbpm]]-model.MUX[model.indx["HAC"]]
#print "model phase x ",model.MUX[model.indx["HAC"]],model.MUX[model.indx[nbpm]]
#phasemodely=model.MUY[model.indx[nbpm]]-model.MUY[model.indx["VAC"]]
#print "model phase x ",model.MUY[model.indx["VAC"]],model.MUY[model.indx[nbpm]]
#drivenphasex=getdrivenphase(phasemodelx*2*pi,Qx,lambdadx,drivetunex,Qx)
#drivenphasey=getdrivenphase(phasemodely*2*pi,Qy,lambdady,drivetuney,Qy)
#print lambdadx,lambdady,Qx,Qy,phasemodelx,phasemodely
#print "Horizontal drive tune ",drivetunex," free oscillation tune ",Qx
#print "Vertical drive tune ",drivetuney," free oscillation tune ",Qy
#print "Horizontal drivenphase ",drivenphasex," and vertical drivenphase ",drivenphasey

htwiss=callandrunmadx(options.accel,"H",options.output,options.bb,Qx,drivetunex,options.madx_var)
data_hor=dosub("H",betax,htwiss)

#sys.exit()

vtwiss=callandrunmadx(options.accel,"V",options.output,options.bb,Qy,drivetuney,options.madx_var)
data_ver=dosub("V",betay,vtwiss)

for bpm in bpms:

    # horizontal
    #phasebpm=phasextot[bpm[1]]
    #print phasebpm
    #print "phasebpm",phasebpm
   # phaseac=phasextot[nbpm]
    #print "total phase ",phaseac,phasebpm
   # drivebeta=betax.BETX[betax.indx[bpm[1]]]
    #drivealfa=betax.ALFX[betax.indx[bpm[1]]]

    #betam=model.BETX[model.indx[bpm[1]]]
    #alfam=model.ALFX[model.indx[bpm[1]]]
    #phasem=model.ALFX[model.indx[bpm[1]]]
    #print betam

    #sac=model.indx["HAC"]
    #s=model.indx[bpm[1]]
    #if(s<sac):
   #     sign="+"
   # else:
    #    sign="-"


        
   # beta=getbeta(lambdadx,phasebpm*2*pi,phaseac*2*pi,drivenphasex,drivetunex,drivebeta,sign)
  #  alfa=getalfa(lambdadx,phasebpm*2*pi,phaseac*2*pi,drivenphasex,drivetunex,drivealfa,sign)
  #  phase=getphase(sign,lambdadx,phasebpm*2*pi,phaseac*2*pi,drivenphasex,drivetunex)



    #print betam

   # horizontal=[bpm, beta,alfa,phase,betam,alfam,phasem,drivenphasex]


   # print >> betaxac, bpm[1], bpm[0],beta, alfa, phase,betam,alfam,phasem,drivenphasex
    print >> betaxac_mo, data_hor[1][bpm[1]][0], data_hor[1][bpm[1]][1], data_hor[1][bpm[1]][2],  data_hor[1][bpm[1]][3],  data_hor[1][bpm[1]][4],0,0,0,0

    #print "x",drivebeta,beta,betam
    
    # vertical
    



   # phasebpm=phaseytot[bpm[1]]
   # phaseac=phaseytot[nbpm]
  #  drivebeta=betay.BETY[betay.indx[bpm[1]]]
  #  drivealfa=betay.ALFY[betay.indx[bpm[1]]]

  #  betam=model.BETY[model.indx[bpm[1]]]
 #   alfam=model.ALFY[model.indx[bpm[1]]]
 #   phasem=model.ALFY[model.indx[bpm[1]]]

  #  sac=model.indx["VAC"]
  #  s=model.indx[bpm[1]]
  #  if(s<sac):
  #      sign="+"
  #  else:
  #      sign="-"
        
  #  beta=getbeta(lambdady,phasebpm*2*pi,phaseac*2*pi,drivenphasey,drivetuney,drivebeta,sign)
  #  alfa=getalfa(lambdady,phasebpm*2*pi,phaseac*2*pi,drivenphasey,drivetuney,drivealfa,sign)
   # phase=getphase(sign,lambdady,phasebpm*2*pi,phaseac*2*pi,drivenphasey,drivetuney)



  #  vertical=[bpm, beta,alfa,phase,betam,alfam,phasem,drivenphasey]

    #print "y", drivebeta,beta,betam

   # print >> betayac, bpm[1], bpm[0],beta, alfa, phase,betam,alfam,phasem,drivenphasey
    print >> betayac_mo,data_ver[1][bpm[1]][0], data_ver[1][bpm[1]][1],data_ver[1][bpm[1]][2], data_ver[1][bpm[1]][3], data_ver[1][bpm[1]][4],0,0,0,0
    
    #print vertical
betaxac.close()
betayac.close()
betaxac_mo.close()
betayac_mo.close()
