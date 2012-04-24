# creating BPMpos3.txt for convert

f=open('install_bpms','r').readlines()
o=open('BPMpos3.txt','w')

text=[]

for line in f:
    if "BPM" in line:
        lineS=line.split()
        lineN=lineS[2],lineS[5]
        text.append(lineN)
        print lineN
        o.writelines(' '.join(lineN))
        o.writelines('\n')
        
        
o.close

 

