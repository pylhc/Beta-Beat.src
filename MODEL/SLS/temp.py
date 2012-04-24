execfile('AllLists.py')
f2=open('quadlocation.txt','w')

var=varSLS()
for i in range(0,len(var)):
    f1=open('temp.txt','r')
    for line in f1:
        varwok=var[i][1:len(var[i])]
        print i,varwok
        if varwok in line:
            sline=line.split()
            f2.write(varwok+' '+sline[3]+'\n')
    f1.close()


f2.close()
            
    
