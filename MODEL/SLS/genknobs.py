fdat=open('SLSring.num.str','r')
fo=open('SLSknobs.py','w')

fo.write('def varSLS():\n')
i=1
for line in fdat:
    sline=line.split()
    try:
        fo.write("\t'"+sline[0]+"',\n")
    except:
        0.0

        
fdat.close()
fo.close()
