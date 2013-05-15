from metaclass import twiss


def iscorr(name):
    if "MDH" in name:
        return 1
    else:
        return 0
    

def checkphaseandbeta(x, k,l):
    #print x.BETX[k]/x.BETX[l],  abs(x.MUX[k]-x.MUX[l])
    if (0.95 < x.BETX[k]/x.BETX[l] < 1.05) and (0.47  < abs(x.MUX[k]-x.MUX[l]) < 0.57 ):
        return 1
    else:
        return 0
    
a=twiss('67ph.opt/twiss.tfs')


pairs=[]
corr1=[]
corr2=[]

corrs=[]



vars={}
for i in range(len(a.NAME)):
    if iscorr(a.NAME[i]):
        corrs.append(i)
        vars[a.NAME[i]]=""
        
vcount=0
for i in range(len(corrs)-3):
    name0=a.NAME[corrs[i]]
    name1=a.NAME[corrs[i+1]]
    name2=a.NAME[corrs[i+2]]
    name3=a.NAME[corrs[i+3]]
    

    if checkphaseandbeta(a, corrs[i], corrs[i+1]):
        pairs.append([corrs[i], corrs[i+1]])
        vars[name0]=vars[name0]+" v"+str(vcount)+" "
        vars[name1]=vars[name1]+" v"+str(vcount)+" "
        vcount+=1
    elif checkphaseandbeta(a, corrs[i], corrs[i+2]):
        pairs.append([corrs[i], corrs[i+2]])
        vars[name0]=vars[name0]+" v"+str(vcount)+" "
        vars[name2]=vars[name2]+" v"+str(vcount)+" "
        vcount+=1
    elif checkphaseandbeta(a, corrs[i], corrs[i+3]):
        pairs.append([corrs[i], corrs[i+3]])
        vars[name0]=vars[name0]+" v"+str(vcount)+" "
        vars[name3]=vars[name3]+" v"+str(vcount)+" "
        vcount+=1
    else:
        print "no match for ", a.NAME[corrs[i]], vcount

for i in range(len(corrs)):
    name=a.NAME[corrs[i]]
    print name, vars[name]
            
        
