import os, datetime, time


now = datetime.datetime.now()
day= now.strftime("%Y-%m-%d")
hour=now.strftime("%H:%M:%S.000")



from optparse import OptionParser
parser = OptionParser()

parser.add_option("-l", "--last",
                help="Last N hours. Use either -l or -t. Priority is given to -l",
                metavar="HOURS", default="1", dest="lastHours")

parser.add_option("-t", "--times",
                help="### NOT yet ready ###  Day and start and end times separated by , ex: 00:00:00.000,01:01:00.000. User either -l or -t",
                metavar="TIMES", default="0", dest="time")

parser.add_option("-s", "--seconds",
                help="Every How many seconds, default 10",
                metavar="SECS", default="60", dest="seconds")


parser.add_option("-o", "--output",
                help="Output path",
                metavar="OUT", default="./", dest="output")



parser.add_option("-c", "--conffile",
                help="Configuration file, normally just use default",
                metavar="CONFFILE", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/TIMBER/ldb.conf", dest="conffile")


parser.add_option("-e", "--executable",
                help="Executable, normally just use default",
                metavar="EXE", default="/afs/cern.ch/group/si/slap/bin/cern-ldb", dest="exe")


(options, args) = parser.parse_args()



#startfrom= now  + datetime.timedelta(hours=-int(options.lastHours))
startfrom= now  + datetime.timedelta(seconds=-60)


startday = startfrom.strftime("%Y-%m-%d")
starthour = startfrom.strftime("%H:%M:%S.000")

print starthour, hour



ACDipole="MKQH.UA43.B1.ACD:FREQUENCY,MKQV.UA43.B1.ACD:FREQUENCY,MKQH.UA43.B2.ACD:FREQUENCY,MKQV.UA43.B2.ACD:FREQUENCY,LHC.BQBBQ.UA47.FFT1_B1:TUNE_H,LHC.BQBBQ.UA47.FFT1_B1:TUNE_V,LHC.BQBBQ.UA43.FFT1_B2:TUNE_H,LHC.BQBBQ.UA43.FFT1_B2:TUNE_V"

#Extract AC dipole information
outfile=options.output+"/ACdipoleWatch.csv"
CommandString=options.exe+" -C "+options.conffile+" -vs \""+ACDipole+"\""+" -t1 \""+startday+" "+starthour+"\" -t2 \""+day+" "+hour+"\" -sa REPEAT -ss "+options.seconds+" -si SECOND   -N "+outfile+ ' '


ch=['QH2','QV2','QH1','QV1','FH1','FH2','FV1','FV2']


while 1:
    time.sleep(5)
    os.system(CommandString)
    fin=open('ACdipoleWatch.csv','r')
    skip=1
    i=0
    for line in fin:
        sline=line.split(',')
        if '-' in line:
            if skip==1:
                skip=0
            else:
                skip=1
                if i<=3:
                    var=ch[i]+'='+sline[1]
                else:
                    s1=sline[1]
                    s1=s1[0:(len(s1)-1)]
                    var=ch[i]+'='+s1+'/11245.0'
                exec(var)
                i=i+1
    #print QH2,QV2,QH1,QV1,FH1,FH2,FV1,FV2
    #sys.exit()

    if abs(QH1-FH1)<0.004:
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 horizontal frequency is dangerous !!!!'
        

    if abs(QH2-FH2)<0.004:
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 horizontal frequency is dangerous !!!!'

        

    if abs(QV1-FV1)<0.004:
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 1 vertical frequency is dangerous !!!!'
        

    if abs(QV2-FV2)<0.004:
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Beam 2 vertical frequency is dangerous !!!!'

        

