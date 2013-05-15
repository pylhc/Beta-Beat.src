from metaclass import twiss


model=twiss('twiss.dat')

data=open('BPMpos3.txt','r')

lines=data.readlines()

names=model.NAME

dictie={}
dictim={}

S=[]

filefile=open("mydictionary.py","w")


print >> filefile, "dictionary={"

for line in lines:

    split=line.split(" ")

    S.append(split[1].replace("\n",""))

    dictie[split[1].replace("\n","")]=split[0]

for name in names:

    if "BP" in name:

        s=model.S[model.indx[name]]
        print s
        dictim[s]=name

print len(dictie),len(dictim)


for ss in S:

    #print ss

    try:
        namem=dictim[float(ss)]
    except:
        print "not in model ",dictie[ss]


    namee=dictie[ss]

    bpm1="SPS."+namem+".H"
    bpm2="SPS."+namem+".V"

    print >> filefile,'    "'+namem+'":["'+bpm1+'",  "'+bpm2+'"] ,'

   

    
print >> filefile, '   "c":"c"}'



filefile.close()
