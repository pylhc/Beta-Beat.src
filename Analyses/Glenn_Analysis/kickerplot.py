###### this python script plots the kick against amp_20,amp_30, Q_x_Y
#
# @ author Glenn Vanbavinckhove
#
#
from metaclass import *
from optparse import *
import os

# option parser
parser = OptionParser()
parser.add_option("-f", "--file base",
                      help="Path to file with twiss file",
                      metavar="FILE", default="./",dest="FILE")

parser.add_option("-m", "--maximum_kick",
                      help="Path to file with twiss file",
                      metavar="MAX", default="./",dest="MAX")
parser.add_option("-b", "--which bpm to look at",
                      help="Path to file with twiss file",
                      metavar="BPM", default="./",dest="BPM")
parser.add_option("-r", "--model",
                      help="Path to file with twiss file",
                      metavar="MODE", default="0",dest="MODE")
(opt, args) = parser.parse_args()

#### command: python kickerplot.py -f cern11 -m 1,2,34 -b BPM001 -r 0

# main

file_base=opt.FILE
model=opt.MODE
print model
max_kick=opt.MAX.split(",")
bpm=opt.BPM
temp=open('./tempplot','w')
temp.write('@ TEXT '+file_base+'\n')
temp.write('* KICK    QX  QXM  QXRMS  QXRMSM   QY  QYM    QYRMS QYRMSM   AMP_30 AMP_30M  AMP_20 AMP_20M \n')
temp.write('$  %le    %le  %le   %le  %le  %le %le    %le   %le     %le %le  %le     %le \n')
for i in range(0,len(max_kick)):


    dir_c='./'+str(file_base)+'_'+str(max_kick[i])+'kV_4kV.txtDir/' # change this if names are different
    file_c=dir_c+str(file_base)+'_'+str(max_kick[i])+'kV_4kV.txt'

    #getting kickdata

    print file_c
    kickdata=twiss(dir_c+'/kickdata')
    linx=twiss(file_c+'_linx')
    liny=twiss(file_c+'_liny')
    mox=twiss(dir_c+'/ALLBPMs_linx')
    moy=twiss(dir_c+'/ALLBPMs_liny')
    qx=linx.Q1
    qxM=mox.Q1
    qy=liny.Q2
    qyM=moy.Q2
    qxrms=linx.Q1RMS
    qxrmsM=mox.Q1RMS
    qyrms=liny.Q2RMS
    qyrmsM=moy.Q2RMS
    amp30=linx.AMP_30[linx.indx[bpm]]
    amp30M=mox.AMP_30[linx.indx[bpm]]
    amp20=linx.AMP_20[linx.indx[bpm]]
    amp20M=mox.AMP_20[linx.indx[bpm]]
    kick2=kickdata.KICK
    kick1=kick2[0]

    if i==0:
        fitqx=qx/kick1
        fitqy=qy/kick1
        fitqamp2=amp20/kick1
        fitqamp3=amp30/kick1
        
        
    temp.write(str(kick1)+'  '+str(qx)+' '+str(qxM)+' '+str(qxrms)+' '+str(qxrmsM)+' '+str(qy)+' '+str(qyM)+' '+str(qyrms)+' '+str(qyrmsM)+' '+str(amp30)+' '+str(amp30M)+' '+str(amp20)+' '+str(amp20M)+'\n')


temp.close()

# gnuplot

# for fitting => linear fit y=a*x+yo  => yo=0 start at 0


GNUcommand=open('plot_kick','w')
GNUcommand.write('set terminal postscript enhanced color solid 20\n')

# for amp30
GNUcommand.write('set xlabel "kick[mm]"\n')
GNUcommand.write('set ylabel "AMP[mm]"\n')
GNUcommand.write('set title "amp_30"\n')
GNUcommand.write('set output "./'+str(file_base)+'_amp30.eps"\n')
GNUcommand.write('l(x)='+str(fitqamp3)+'*x\n')
if model=='1':
    GNUcommand.write('p [][]"./tempplot" u 1:10 t \'real\' w l, "" u 1:11 t \'simul\' w l\n') # with model
elif model=='0':
    GNUcommand.write('p [][]"./tempplot" u 1:10 t \'real\',l(x) \n')

# for amp20
GNUcommand.write('set xlabel "kick [mm]"\n')
GNUcommand.write('set ylabel "AMP[mm]"\n')
GNUcommand.write('set title "amp_20"\n')


GNUcommand.write('set output "./'+str(file_base)+'_amp20.eps"\n')

if model=='1':
    GNUcommand.write('p [:][0:15]"./tempplot" u 1:12 t \'real\', "" u 1:13 t \'simul\'\n')
elif model=='0':
    GNUcommand.write('f(x) = a*x**2 + b*x + c\n')
    GNUcommand.write('p [][]"./tempplot" u 1:12 t \'real\'\n')
    GNUcommand.write('fit f(x) "./tempplot" using 1:12 via a,b,c\n')
    GNUcommand.write('replot')
# for Qx
GNUcommand.write('set xlabel "kick [mm]"\n')
GNUcommand.write('set ylabel "Qx[-]"\n')
GNUcommand.write('set title "amp_20"\n')


GNUcommand.write('set output "./'+str(file_base)+'_Qx.eps"\n')
GNUcommand.write('l(x)='+str(fitqx)+'*x\n')
if model=='1': 
    GNUcommand.write('p [][]"./tempplot" u 1:2:4 w errorbars t \'real\', "" u 1:3:5 w errorbars t \'simul\'\n')
elif model=='0':
    GNUcommand.write('p [][]"./tempplot" u 1:2:4 w errorbars t \'real\',l(x)\n')

# for Qy
GNUcommand.write('set xlabel "kick [mm]"\n')
GNUcommand.write('set ylabel "Qy[-]"\n')
GNUcommand.write('set title "amp_20"\n')


GNUcommand.write('set output "./'+str(file_base)+'_Qy.eps"\n')
GNUcommand.write('l(x)='+str(fitqx)+'*x\n')
if model=='1':
    GNUcommand.write('p [][]"./tempplot" u 1:6:8 w errorbars t \'real\', "" u 1:7:9 w errorbars t \'simul\'\n')
elif model=='0':
    GNUcommand.write('p [][]"./tempplot" u 1:6:8 w errorbars t \'real\',l(x) \n')
    

GNUcommand.close()

os.system('gnuplot  plot_kick')
