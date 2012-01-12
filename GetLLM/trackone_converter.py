
f=open('trackone','r')
g=open('trackone_c','w')

elemno=529 # enter no of bpm +1 because of "end marker"
turnno=256

i=1
lsi=[]
ldi=[]
ld=[]
te=[]
for line in f:
    if i<=55:
        0.0
    elif i%2==0:
        lsi=line.split()
        bpm=lsi[5]
    elif i%2==1:
        ldi=line.split()
        ld.append([bpm,ldi[2],ldi[4],ldi[8]])
    if i>55 and i!=56 and (i-55)%(elemno*2)==0:
        te.append(ld)
        ls=[]
        ld=[]
    i+=1


f.close

g.write('title \n')

for elem in range(0,elemno-1):
    for xy in range(0,2):
        g.write(str(xy)+' '+te[0][elemno-2-elem][0]+' '+te[0][elemno-2-elem][3]+' ')
        for turn in range(0,turnno):
            #print elem,turn,xy
            g.write(str(te[turnno-1-turn][elemno-2-elem][xy+1]+' '))
        g.write('\n')

