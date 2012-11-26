import json
import os,sys


def convert(fin,fout):
    print "Converting",fin,"-->",fout
    ret={}
    mod=fin[:-3]
    exec('import '+mod)
    exec('a=dir('+mod+')')
    for key in a:
        if '_' not in key:
            print '  ',key
            exec('ret[key]='+mod+'.'+key+'()')
    file(fout,'w').write(json.dumps(ret,sort_keys=True,indent=2))

for f in os.listdir('.'):
    if f[:8]=='AllLists' and f[-3:]=='.py':
        fout=f[:-2]+'json'
        fin=f
        convert(fin,fout)


