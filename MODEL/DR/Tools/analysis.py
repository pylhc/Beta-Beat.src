from metaclass import twiss

# load model twiss 
model=twiss('dr.tfs')

MUX=model.MUX
MUY=model.MUY

names=model.NAME
S=model.S

data=open('phases.tfs','w')

# for BPMs
#bpmdef=open('bpm_def.madx','w')
#bpminstall=open('bpm_install.madx','w')

print >> data, "NAME NAME2 S S2 PHX PHY"
print >> data, "%s %s %le %le %le %le"

for i in range(len(names)-1):

    name1=names[i]
    name2=names[i+1]
    loc=S[i]
    loc2=S[i+1]
    phx=(MUX[i+1]-MUX[i])%1
    phy=(MUY[i+1]-MUY[i])%1

    print >> data,name1,name2,loc,loc2,phx,phy

#for i in range(len(names)):

#    el="bpm"+str(i)
#    pos=S[i]

#    print >> bpmdef,el," : monitor , l=0 ;"
 #   print >> bpminstall,"install,element=",el,"at = ",pos,";"


#bpmdef.close()
#bpminstall.close()
data.close()
