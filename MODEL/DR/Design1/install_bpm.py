from metaclass import twiss
import sys

model=twiss('quad.tfs')

bpminstall=open('bpms_install.madx','w')
bpmseq=open('bpms_seqedit.madx','w')
bpmptc=open('bpms_ptc.madx','w')

names=model.NAME
stot=model.S


# phase advances
ph=[] # pos
phx=0
phy=0
summ=0

for i in range(len(names)-1):

    name1=names[i]
    name2=names[i+1]
    s=stot[i]

    phx=(model.MUX[i+1]-model.MUX[i]+phx)%1
    phy=(model.MUY[i+1]-model.MUY[i]+phy)%1

    #print model.MUX[i+1]-model.MUX[i],model.NAME[i],model.NAME[i+1],model.S[i],model.S[i+1]

    #print "phx ",phx

    if (phx<0.12) or (phx>0.35): # 45 and 135

        #print phx
        
        keyfi=i
        phx=phx
        phy=phy

        summ=1
        passs=0

     #   print "small", phx

    else:

        #print "pass",phx

        if summ==1:
            keyfirst=keyfi

        else:
            keyfirst=i
        
        keylast=i+1
        passs=1
        phx=0
        phy=0
        keyfi=0
        summ=0


    if passs==1:

        ph.append(keyfirst)

        print "adding",(model.MUX[keylast]-model.MUX[keyfirst])%1

    #print phx,phy

#sys.exit()

count=0

for key in ph:

    s=model.S[key]

    print >> bpmseq, "install,element= bpm"+str(count)+" ,at =  "+str(s)+" ;"
    print >> bpminstall,"bpm"+str(count)+" : monitor , l="+ str(0)+" ;"
    print >> bpmptc,"ptc_observe, place=\"bpm"+str(count)+"\";"

    count=count+1
    
bpminstall.close()
bpmseq.close()
bpmptc.close()
