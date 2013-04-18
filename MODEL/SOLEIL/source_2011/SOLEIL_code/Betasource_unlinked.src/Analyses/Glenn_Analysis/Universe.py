###########################################
#
# Created by Glenn Vanbavinckhove
# Date : 01/04/09
#
###########################################
#
#
# Version 0.1 :
#               => how to run : python universe.py -b ./Beta-beat.src  -f file2Data -o ./  s- 0

# import
from metaclass import twiss
from optparse import *
import os,sys
import time
import shutil
import math



# option parser
parser = OptionParser()
parser.add_option("-b", "--source", 
                      help="Path to the beta-beat source directory",
                      metavar="SOURCE", default="./",dest="SOURCE")
parser.add_option("-f", "--file",
                      help="Path to file with twiss file",
                      metavar="PATH", default="./",dest="PATH")
parser.add_option("-o", "--outputfile",
                      help="Where to put the stuff",
                      metavar="OUTPUT", default="./",dest="OUTPUT")
parser.add_option("-m", "--MAD",
                      help="path to MAD FILES",
                      metavar="MAD", default="0",dest="MAD")
    
(opt, args) = parser.parse_args()

######### Main part
sourceDIR=opt.SOURCE
outputDir=opt.OUTPUT
MADP=opt.MAD
twissFILE=twiss(str(opt.PATH))

dirp=twissFILE.DIR
files=twissFILE.FILE
TunesX=twissFILE.QX
TunesY=twissFILE.QY
Kick=twissFILE.KICK
Turns=twissFILE.TURNS
HKick=twissFILE.HK
VKick=twissFILE.VK
Keyword=twissFILE.KEYWORD
Relation=twissFILE.REL
filess=[]

print 'pre-check of files will be carried out ... '
print 'WARNING: if file does not exist system will exit !!'
for i in range(0,len(files)):
    sourcefile=dirp+'/'+files[i]
    check=os.path.exists(sourcefile)

    if check==False:
        print str(sourcefile)+' file does not exist... system exit'
        sys.exit()

print 'good... all files exists'

print 'looking for highest relation'
xmax=0

#dell=[]
#filenum=[]
#number=[]
#for i in range(0,len(Relation)):

    
    
    #filterr=filter(lambda x: x == Relation[i], Relation)
      
    #if filterr[0]==Relation[i]:
           # number.append(i)
    #filenum.append(number)
    
    #if Relation[i]!=0: 
        #try:
            #dell.index(filterr)
        #except:
            #dell.append(filterr)
 
#print number
#print filenum
#sys.exit()

#for i in range(0,len(Relation)):
    
    #sourcefile=dirp+'/'+files[i]
    #if Relation[i]>xmax:
        #xmax=Relation[i]
   

#print 'Number of relations found '+str(xmax)

#rela=[]
#numrel=[]

hippie=['d']

for i in range(0,len(hippie)):

    # copying file into dir
    sourcefile=dirp+'/'+files[i]
    destination=outputDir+files[i]+'Dir'
    destinationfile=destination+'/'+files[i]
    os.mkdir(destination)
    shutil.copyfile(sourcefile,destinationfile)
    filess.append(destinationfile.replace)


    #drive.inp
    driveFile=open(destination+'/Drive.inp','w')
    driveFile.write('KICK='+str(Kick[i])+'\n')
    driveFile.write('CASE(1[H], 0[V])=1\n')
    driveFile.write('KPER(KICK PERCE.)=0.5\n')
    driveFile.write('TUNE X ='+str(TunesX[i])+'\n')
    driveFile.write('TUNE Y ='+str(TunesY[i])+'\n')
    driveFile.write('PICKUP START=0\n')
    driveFile.write('PICKUP END=550\n')
    driveFile.write('NORMALISATION(1[yes])=1\n')
    driveFile.write('ISTUN =0.02\n')
    driveFile.write('BETABEATING(1[yes])=0\n')
    driveFile.write('IR(1[Real], 0[Comp])=0\n')
    driveFile.write('LABEL RUN (1[yes])=0\n')
    driveFile.write('FORMAT (0[SPS],1[HERA],2[RHIC])=2\n')
    driveFile.write('NOISEPATH =noisefiles/\n')
    driveFile.write('WINDOWa1=0.01\n')
    driveFile.write('WINDOWa2=0.07\n')
    driveFile.write('WINDOWb1=0.35\n')
    driveFile.write('WINDOWb2=0.45\n')

    driveFile.close()

    #drivingterms
    termsFile=open(destination+'/DrivingTerms'+'','w')
    termsFile.write(str(destinationfile)+'  0   '+str(Turns[i]))
    termsFile.close()



    #Run Drive God lin
    os.system(sourceDIR+'/DRIVE_src/Drive_God_lin   '+destination+'/')



    #Run SVD
    os.system(sourceDIR+'//SVD/phSVD32bits/phSVD.py   -p '+destination+'/'+' -a SOLEIL -m 0 -f '+files[i])


    #Run gnuplot script for frequency and TBT
    gnuS=destination+'/plot4Universe'
    GNUcommand=open(gnuS,'w')
    GNUcommand.write('set terminal postscript enhanced color solid 20\n')
    GNUcommand.write('set xlabel "s [m]"\n')
    GNUcommand.write('set ylabel "position[mm]"\n')
    GNUcommand.write('set title "'+str(Keyword[i])+'"\n')
    GNUcommand.write('set logscale y\n')
    GNUcommand.write('unset key\n')
    GNUcommand.write('set output "'+str(destination)+'/freqx.eps"\n')

    # freq x
    GNUcommand.write('p [-0.005:][]"'+destination+'/bpm001.x" u 1:2 w i\n') # change name here for accelerator

    
    # freq y
    GNUcommand.write('set output "'+str(destination)+'/freqy.eps"\n')
    GNUcommand.write('p [-0.005:][]"'+destination+'/bpm001.y" u 1:2 w i\n') # change name here for accelerator
    GNUcommand.write('\n')

    GNUcommand.write('unset logscale y\n')
    GNUcommand.write('\n')

    # TBT x
    os.system('grep \'0  bpm001\' '+destinationfile+' | awk \'BEGIN{RS=" "} ($1*1.0)**2>0{print}\'  > tx\n')
    GNUcommand.write('set output "'+str(destination)+'/TBTx.eps"\n')
    GNUcommand.write('p "tx"\n') # change name here for accelerator
    GNUcommand.write('\n')

    # TBT y
    os.system('grep \'1  bpm001\' '+destinationfile+' | awk \'BEGIN{RS=" "} ($1*1.0)**2>0{print}\'  > ty\n')
    GNUcommand.write('set output "'+str(destination)+'/TBTy.eps"\n')
    GNUcommand.write('p "ty"\n') # change name here for accelerator
    
    
    GNUcommand.write('\n')

    # gnuplot for noise
    GNUcommand.write('set output "'+str(destination)+'/noiseSVD.eps"\n')
    GNUcommand.write('set size 1.05,0.95\n')
    GNUcommand.write('set pointsize 1\n')
    GNUcommand.write('set key top samplen 1\n')
    GNUcommand.write('\n')

    GNUcommand.write('set ylabel "Counts" 0.8\n')
    GNUcommand.write('set xlabel "BPM rms noise [mm]"\n')
    GNUcommand.write('set ytics 20\n')
    GNUcommand.write('bx=0.14 #--- approx scaling factor for mm ---> deg\n')
    GNUcommand.write('by=0.17 #--- approx scaling factor for mm ---> deg\n')
    GNUcommand.write('p "'+destinationfile+'_nhistx" u ($1/1e3):($2/2.38) t "Horizontal" w histeps lt 1 lw 4 ,\\\n')
   
    GNUcommand.write('"'+destinationfile+'_nhisty" u ($1/1e3):($2/1.96) t "Vertical" w histeps lt 3 lw 4 \n')

    GNUcommand.close()

    os.system('gnuplot  '+gnuS)

    # kick calculation
    kickF=destination+'/kickdata'
    kicker=open(kickF,'w')
    
    betax=10.87
    betaxbpm=13.83
    betay=7.99
    betaybpm=12.04

    kicker.write('* PLANE         BETAS        BETABPM      FACTOR    VOLTAGE      KICK \n')
    kicker.write('$  %s            %le            %le        %le        %le         %le \n')

    factorx=1.7
    factory=0.58
    xx=math.sqrt(betax/betaxbpm)*factorx*HKick[i]/1000

    yy=math.sqrt(betay/betaybpm)*factory*VKick[i]/1000

    kicker.write('H         '+str(betax)+'        '+str(betaxbpm)+'      '+str(factorx)+'    '+str(HKick[i])+'      '+str(xx)+'\n')
    kicker.write('V         '+str(betay)+'        '+str(betaybpm)+'      '+str(factory)+'    '+str(VKick[i])+'      '+str(yy)+'\n')
    kicker.close()

    # commad for madx
    #shutil.copyfile(MADP+'/AddBpmErrors',destination+'/AddBpmErrors')
    #os.system('chmod 777  '+destination+'/AddBpmErrors')
    madxS=destination+'/job.'+files[i]+'.madx'
    madxFile=open(madxS,'w')
    madxFile.write('BEAM ,PARTICLE=electron,ENERGY=0.27505112E+01;\n')
    madxFile.write('call, file="'+MADP+'/sol.2009.madx";\n')
    madxFile.write('USE , sequence=ring;\n')
    madxFile.write('call, file="'+MADP+'/define_bpms_new";\n')
    madxFile.write('seqedit, sequence=ring;\n')
    madxFile.write('call, file="'+MADP+'/install_bpms_new";\n')
    madxFile.write('endedit;\n')
    madxFile.write('\n')
    madxFile.write('use, sequence=ring;\n')
    madxFile.write('\n')
    madxFile.write('!!!!!!!!!!!!!!!!!!!!!!!!\n')
    madxFile.write('! Match tunes\n')
    madxFile.write('!!!!!!!!!!!!!!!!!!!!!!!!\n')
    madxFile.write('\n')
    madxFile.write('match;\n')
    madxFile.write('vary, name=kqp7, step=0.00000000001;\n')
    madxFile.write('vary, name=kqp9, step=0.00000000001;\n')
    madxFile.write('constraint, range=#e, mux='+str(TunesX[i])+', muy='+str(TunesY[i])+';\n')
    madxFile.write('lmdif, tolerance=1e-6;\n')
    madxFile.write('endmatch\n')
    madxFile.write('\n')
    madxFile.write('! Errors in SHORT quads\n')
    madxFile.write('SELECT,FLAG=ERROR, clear;SELECT,FLAG=ERROR, PATTERN="QP1";SELECT,FLAG=ERROR, PATTERN="QP3";SELECT,FLAG=ERROR, PATTERN="QP4";SELECT,FLAG=ERROR, PATTERN="QP5";SELECT,FLAG=ERROR, PATTERN="QP6";SELECT,FLAG=ERROR, PATTERN="QP8";\n')
    madxFile.write('EFCOMP, ORDER:=1,RADIUS:=0.03,DKNR:={0,0,-1.6e-4 , -3.4e-4 ,0 , 2.4e-4},DKSR:={0,0,0.5e-4};\n')
    madxFile.write('\n')
    madxFile.write('! Errors in LONG quads\n')
    madxFile.write('SELECT,FLAG=ERROR, clear;SELECT,FLAG=ERROR, PATTERN="QP2";SELECT,FLAG=ERROR, PATTERN="QP7";\n')
    madxFile.write('EFCOMP, ORDER:=1,RADIUS:=0.03,DKNR:={0,0,2.9e-4 , -8.6e-4 ,0 , 0.7e-4},DKSR:={0,0,0.5e-4};\n')
    madxFile.write('\n')
    madxFile.write('SELECT,FLAG=ERROR, class=quadrupole;\n')
    madxFile.write('esave, file=err;\n')

    madxFile.write('PTC_CREATE_UNIVERSE;PTC_CREATE_LAYOUT, model=2, method=6, nst=10;\n')
    madxFile.write('call, file="'+MADP+'/observe_bpms_new";\n')
    madxFile.write('PTC_START, x='+str(xx)+', y='+str(yy)+';\n')
    madxFile.write('PTC_TRACK, deltap=0.0000, icase=5, turns=2000,dump, onetable;\n')
    madxFile.write('PTC_TRACK_END;\n')
    madxFile.write('PTC_END;\n')
    madxFile.write('\n')

    #madxFile.write('system, "cp ./trackone  '+destination+'/trackone"\n')
    #madxFile.write('system, "AddBpmErrors"\n')

    madxFile.write('\n')
    madxFile.write('stop;\n')

    madxFile.close()

    print 'madxp < '+str(madxS)
    os.system('madxp <  '+madxS)

    #shutil.copyfile('./job.'+files[i]+'.madx',destination+'/job.'+files[i]+'.madx')
    #shutil.copyfile('./trackone',destination+'/trackone')
    #os.system('chmod 777  '+destination+'/trackone')
    #os.remove('./job.'+files[i]+'.madx')
    #os.remove('./trackone')
    os.system('AddBpmErrors')
    shutil.copyfile('./ALLBPMs',destination+'/ALLBPMs')
    
    #drivingterms
    os.remove(destination+'/DrivingTerms')
    termsFile=open(destination+'/DrivingTerms'+'','w')
    termsFile.write(str(destination)+'/ALLBPMs'+'  0   '+str(Turns[i]))
    termsFile.close()

    # running drive god_lin second time
    os.system(sourceDIR+'/DRIVE_src/Drive_God_lin   '+destination+'/')



    # running gnuplot for simul and real
    realx=destinationfile+'_linx'
    #amprx=realx.AMP_20
    simulx=destination+'/ALLBPMs_linx'
    #ampsx=simulx.AMP_20
    #realy=twiss(destinationfile+'_liny')
    #simuly=twiss(destination+'ALLBPMs_liny')
    
    gnuS2=destination+'/plot4Universe2'
    GNUcommand2=open(gnuS2,'w')
    GNUcommand2.write('set terminal postscript enhanced color solid 20\n')
    GNUcommand2.write('set xlabel "s [m]"\n')
    GNUcommand2.write('set ylabel "Amplitude20/Amplitude10[mm]"\n')
    GNUcommand2.write('set title "'+str(Keyword[i])+'"\n')

    GNUcommand2.write('set output "'+str(destination)+'/amp_sim_re_x.eps"\n')

    GNUcommand2.write('p "'+realx+'" u 2:($14*6) t \'real\' w l, "'+simulx+'" u 2:14 t \'simul\' w l')
    
    GNUcommand2.close()

    os.system('gnuplot '+gnuS2)

    d=i+1
    print str(d)+' is finished'
