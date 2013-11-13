#
# => getting fit data, replacing all the stufffff
#
#
#

from metaclass import *
import os
from math import cos,sin
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-s", "--confile",
                help="file that defines segment",
                metavar="config", default="none",dest="config")
parser.add_option("-p", "--path",
                help="path to measurement",
                metavar="path", default="./",dest="path")
parser.add_option("-d", "--store",
                help="where to store dir",
                metavar="dir", default="./",dest="dir")
parser.add_option("-t", "--tune",
                help="tune for fit ",
                metavar="tune", default="5",dest="tune")
parser.add_option("-a", "--amp",
                help="amp guess for hor and ver",
                metavar="amp", default="5",dest="amp")
parser.add_option("-f", "--phase",
                help="phase (fase in dutch) guess for hor and ver",
                metavar="pha", default="5",dest="pha")
parser.add_option("-r", "--repository",
                help="beta beat location",
                metavar="bb", default="./",dest="bb")
parser.add_option("-c", "--cut",
                help="cut",
                metavar="cut", default="0.05",dest="cut")
parser.add_option("-g", "--gnuplot",
                help="gnuplot",
                metavar="gnuplot", default="gnuplot",dest="gnuplot")
parser.add_option("-m", "--plane",
                help="transverse plane",
                metavar="plane", default="H",dest="plane")
parser.add_option("-o", "--option",
                help="rerun only function (1) or also fit(0)",
                metavar="runoption", default="0",dest="runoption")


(options, args) = parser.parse_args()


def getdata(path):


    data=open(path+'/log','r')
    lines=data.readlines()
    amp=0;erra=0;phase=0;errp=0
    for line in lines:

        if "a" in line and "+/-" in line:

            linesplit=line.split()
            amp=linesplit[2]
            erra=linesplit[4]
            
        if "p" in line and "+/-" in line:
            linesplit=line.split()
            phase=linesplit[2]
            errp=linesplit[4]

  
    data.close

    #print phase,errp

    return [amp,erra,phase,errp]

# writing gnuplot scripts
def forgnuplot(path,tune,amp,phase,files,name,cpath,cut,dir2store,dir):

    filename=path+'/var4plot.sh'
    file4nad=open(filename,'w')
    #file4nad.write('#!bash \n')
    file4nad.write('sed -e \'s/%TUNE/\''+str(tune)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%AMP/\''+str(amp)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%PHASE/\''+str(phase)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%INPUT/\'\"'+str(files.replace('/','\/'))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%OUTPUT/\''+str(name)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%LOG/\'\"'+str(dir2store.replace('/','\/'))+'\"\'/g\' \\\n')
    file4nad.write('    -e \'s/%CUT/\''+str(cut)+'\'/g\' \\\n')
    file4nad.write('    -e \'s/%DIR/\'\"'+dir+'\"\'/g\' \\\n')
    file4nad.write('<'+cpath+'/PhaseFit/gplot_fit.mask > '+path+'/gplot \n')
    file4nad.close()
    os.system("chmod 777 "+str(filename))
    os.system(str(filename))
    os.system("chmod 777 "+str(path+'/gplot'))
    os.system(options.gnuplot+" "+path+'/gplot')


# filtering bpms to get the correct values
def filterbpms(measuredx,path,startBPM,endBPM,dir):

    if dir=="H":
        filex=open(path+"/phasefitx.out","w")
        filex.write("* NAME  S  MUXMDL  PHASEX  STDPHX  PHXMDL\n")
    else:
        filex=open(path+"/phasefity.out","w")
        filex.write("* NAME  S  MUYMDL  PHASEY  STDPHY  PHYMDL\n")
    names=measuredx.NAME
    in1=measuredx.indx[startBPM]
    in2=measuredx.indx[endBPM]
    
    filex.write("$ %s    %le  %le    %le   %le   %le\n")

    for bpm in names:

        ind=measuredx.indx[bpm]
        
        if  ind > in1 and ind < in2:
            if dir=="H":
                phase=measuredx.PHASEX[ind]
                phasemdl=measuredx.PHXMDL[ind]
                muxmdl=measuredx.MUXMDL[ind]
                stdph=measuredx.STDPHX[ind]
                S=measuredx.S[ind]
            else:
                phase=measuredx.PHASEY[ind]
                phasemdl=measuredx.PHYMDL[ind]
                muxmdl=measuredx.MUYMDL[ind]
                stdph=measuredx.STDPHY[ind]
                S=measuredx.S[ind]

            filex.write(bpm+" "+str(S)+" "+str(muxmdl)+" "+str(phase)+" "+str(stdph)+" "+str(phasemdl)+"\n")

    filex.close()
                



file4seg=twiss(options.config)
dir2store=options.dir
startBPM=file4seg.SBPM
endBPM=file4seg.EBPM
print "Will perform fit from "+startBPM+" to "+endBPM
bpms=file4seg.NAME
path=options.path
tune=options.tune
print tune
amp=options.amp
phase=options.pha
cut=options.cut
cpath=options.bb
name=file4seg.LABEL
    

if options.plane=="H":
    files=dir2store+"/phasefitx.out"
    measured=twiss(options.path+'/getphasex.out')
    
else:
    files=dir2store+"/phasefity.out"
    measured=twiss(options.path+'/getphasey.out')
filterbpms(measured,dir2store,startBPM,endBPM,options.plane)
infile= open(dir2store+'/getfit_'+name+'_'+options.plane+'.out','w')
names=measured.NAME
if options.runoption=="0":

    print "am performing fit"
    if os.path.exists(dir2store+"/log"):
        os.system("rm "+dir2store+"/log")
    forgnuplot(dir2store,tune,amp,phase,files,name,cpath,cut,dir2store,options.plane)
    fromfit=getdata(dir2store)
    infile.write("@ AMP %le "+str(fromfit[0])+"\n")
    infile.write("@ ERRAMP %le "+str(fromfit[1])+"\n")
    infile.write("@ PHASE %le "+str(fromfit[2])+"\n")
    infile.write("@ ERRPHASE %le "+str(fromfit[3])+"\n")
    
    A=float(fromfit[0])
    P=float(fromfit[2])
    Q=tune
    
else:

    print "am only running function"
    fromfit=getdata(dir2store)
    A=float(amp)
    P=float(phase)
    Q=tune
    infile.write("@ AMP %le "+str(A)+"\n")
    infile.write("@ ERRAMP %le "+str(fromfit[1])+"\n")
    infile.write("@ PHASE %le "+str(P)+"\n")
    infile.write("@ ERRPHASE %le "+str(fromfit[3])+"\n")






infile.write("* NAME  S  MUXMDL  BEAT  STDPH  FIT\n")
infile.write("$ %s    %le  %le    %le   %le   %le\n")

in1=measured.indx[startBPM]
in2=measured.indx[endBPM]

for bpm in names:

    if options.plane=="H":
        beat=measured.PHASEX[measured.indx[bpm]]-measured.PHXMDL[measured.indx[bpm]]
    else:
        beat=measured.PHASEY[measured.indx[bpm]]-measured.PHYMDL[measured.indx[bpm]]
        



    ind=measured.indx[bpm]
        
    if  ind > in1 and ind < in2:

        if beat**2<float(cut)**2:

            if options.plane=="H":
                x=measured.PHXMDL[measured.indx[bpm]]
                y=measured.MUXMDL[measured.indx[bpm]]
                err=measured.STDPHX[measured.indx[bpm]]
                s=measured.S[measured.indx[bpm]]
            else:
                x=measured.PHYMDL[measured.indx[bpm]]
                y=measured.MUYMDL[measured.indx[bpm]]
                err=measured.STDPHY[measured.indx[bpm]]
                s=measured.S[measured.indx[bpm]]
    
            point=A*(sin(2*abs(y-P)*2*pi-2*pi*float(Q))-sin(2*abs(y+x-P)*2*pi-2*pi*float(Q)))
            #print float(Q)
            infile.write(bpm+" "+str(s)+" "+str(y)+" "+str(beat)+" "+str(err)+" "+str(point)+"\n")
        
        else:

            print "skipping "+bpm


infile.close()


print "fit class finished"
