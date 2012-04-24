#! /usr/bin/env python

import Tkinter as Tk
import os, datetime, time



class Frame(Tk.Frame):
    def __init__(self, master=None):
        Tk.Frame.__init__(self, master)
        self.master.title('AC Dipole Watcher')
        self.b1h_col='green'
        self.b1v_col='green'
        self.b2h_col='green'
        self.b2v_col='green'

        self.la=Tk.Label(self, text='Beam1 H', bg=self.b1h_col, font=('Times','100'))
        self.lb=Tk.Label(self, text='Beam1 V', bg=self.b1v_col, font=('Times','100'))
        self.lc=Tk.Label(self, text='Beam2 H', bg=self.b2h_col, font=('Times','100'))
        self.ld=Tk.Label(self, text='Beam2 V', bg=self.b2v_col, font=('Times','100'))
        self.la.pack()
        self.lb.pack()
        self.lc.pack()
        self.ld.pack()

        self.bind_all('<1>',self.starting)

    def starting(self,event):
        self.after(1,self.watching)


    def watching(self):
        time.sleep(1)

        now = datetime.datetime.now()
        day= now.strftime("%Y-%m-%d")
        hour=now.strftime("%H:%M:%S.000")

        startfrom= now  + datetime.timedelta(seconds=-60)


        startday = startfrom.strftime("%Y-%m-%d")
        starthour = startfrom.strftime("%H:%M:%S.000")
        ACDipole="MKQH.UA43.B1.ACD:FREQUENCY,MKQV.UA43.B1.ACD:FREQUENCY,MKQH.UA43.B2.ACD:FREQUENCY,MKQV.UA43.B2.ACD:FREQUENCY,LHC.BQBBQ.UA47.FFT1_B1:TUNE_H,LHC.BQBBQ.UA47.FFT1_B1:TUNE_V,LHC.BQBBQ.UA43.FFT1_B2:TUNE_H,LHC.BQBBQ.UA43.FFT1_B2:TUNE_V"


        outfile="./ACdipoleWatch.csv"
        CommandString="/afs/cern.ch/group/si/slap/bin/cern-ldb -C /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/TIMBER/ldb.conf -vs \""+ACDipole+"\""+" -t1 \""+startday+" "+starthour+"\" -t2 \""+day+" "+hour+"\" -sa REPEAT -ss 60 -si SECOND   -N "+outfile+ ' '



        ch=['QH2','QV2','QH1','QV1','FH1','FH2','FV1','FV2']

        os.system(CommandString)
        fin=open('./ACdipoleWatch.csv','r')
        i=0
        skip=1
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


        if abs(QH1-FH1)>0.004:
            self.b1h_col='green'
        else:
            self.b1h_col='red'

        if abs(QV1-FV1)>0.004:
            self.b1v_col='green'
        else:
            self.b1v_col='red'

        if abs(QH2-FH2)>0.004:
            self.b2h_col='green'
        else:
            self.b2h_col='red'

        if abs(QV2-FV2)>0.004:
            self.b2v_col='green'
        else:
            self.b2v_col='red'

    
        self.la.configure(bg=self.b1h_col)
        self.lb.configure(bg=self.b1v_col)
        self.lc.configure(bg=self.b2h_col)
        self.ld.configure(bg=self.b2v_col)
        self.la.pack()
        self.lb.pack()
        self.lc.pack()
        self.ld.pack()
       
        self.after(100,self.watching)



if __name__=='__main__':
    root= Frame()
    root.pack()    
    root.mainloop()



