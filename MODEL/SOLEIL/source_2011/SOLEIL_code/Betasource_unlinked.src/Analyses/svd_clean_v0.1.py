#--- Newer version of svd_clean with numpy
#--- use rhicdata class and multi-threading
#--- R. Calaga Mar 4, 2011
import sys, time,optparse,os
try:
    import numpy as npy;
    from rhicdata25 import rhicdata
    from handythread import foreach
except ImportError: print "Import Error"; sys.exit()    

#---------------------------------------
def parsse():
    parser = optparse.OptionParser()
    parser.add_option("-f", "--file",
                      help="Specify the file name",
                      metavar="FILE", default="None",dest="FILE")
    parser.add_option("-p", "--peak2peak",
                      help="min peak to peak amplitude [default=0.2mm]",
                      metavar="PK2PK", default=0.2,dest="PK2PK")
    parser.add_option("-s", "--sumsquare",
                      help="amplitude for svd cut [default=0.925]",
                      metavar="SSQ", default=0.925,dest="SSQ")
    parser.add_option("-v", "--svals",
                      help="# of singular vals to retain [default=500]",
                      metavar="SVALS", default=1000,dest="SVALS")
    parser.add_option("-t", "--turns",
                      help="Turn number to start [1...N Turns]",
                      metavar="NT", default=0,dest="NT")
    (opt, args) = parser.parse_args()
    return opt, args

def rBPM(data):
    NT,NB=npy.shape(data)
    if NB>10:
        datavg=npy.average(data,axis=0)
        U,S,V=npy.linalg.svd((data-datavg)/npy.sqrt(NT),\
                             full_matrices=0);
        V=npy.transpose(V);row,col=npy.shape(V);
        nr=[];bidx=[];gab=[];VBD=npy.array([])
        #-- reverse sort
        VS=npy.sort(npy.abs(V),axis=0);VS=VS[::-1]
        VAS=npy.argsort(npy.abs(V),axis=0);VAS=VAS[::-1]
        for k in range(col):
            for j in range(row):
                if npy.add.reduce(VS[:j,k]**2)>float(opt.SSQ):
                    if j<5: bidx.append(k)
                    break
        #-- remove redundancy
        try: VBD=npy.unique1d(npy.take(VAS[0,:],(bidx),1))
        except: pass
        for ee in range(row):
            if ee not in VBD: gab.append(ee)
        return gab
    else: print "Less than 10 bpms, exiting.."; sys.exit()

def clean(data):
    #-- no sqrt(N) norm like old svd_clean
    datavg=npy.average(data,axis=0); norm=npy.sqrt(len(data))
    U,S,V=npy.linalg.svd(data-datavg,full_matrices=0)
    Sn=npy.copy(S); Sn[int(opt.SVALS):]=0.0
    BN=npy.dot(U,npy.dot(npy.diag(Sn),V))
    return BN+datavg

def red(fnm):
    if os.path.isfile(fnm):a=rhicdata(fnm)
    else: print fnm, "doesn't exists"; sys.exit()
    print "Selected file:", opt.FILE    
    tx=npy.array([j.data[int(opt.NT)-1:] for j in a.H])
    ty=npy.array([j.data[int(opt.NT)-1:] for j in a.V])
    return a, npy.transpose(tx),npy.transpose(ty)

def svdClean():
    #-- read file
    t0=time.time(); a,tx,ty=red(opt.FILE)
    #print "File read in",round(time.time()-t0,1),'s'    
    ntx,nbx=npy.shape(tx);nty,nby=npy.shape(ty)
    print '[H, V] bpms: [',nbx, nby,']'
    print '[H, V] turns: [',ntx, nty,']'

    #--- peak-2-peak cut, for LHC convert to microns
    print "Peak to peak cut:",opt.PK2PK, "mm"        
    pkx=npy.nonzero(npy.ptp(tx,axis=0)>float(opt.PK2PK))[0]
    pky=npy.nonzero(npy.ptp(ty,axis=0)>float(opt.PK2PK))[0]
    tx=npy.take(tx,pkx,1);ty=npy.take(ty,pky,1)
    print '[H,V] BPMs after P2P cut:',len(pkx),len(pky)
    
    #--- svd cut
    #t0=time.time()
    #gdx,gdy=foreach(rBPM,[tx,ty],threads=2,return_=True)
    gdx=rBPM(tx);gdy=rBPM(ty); #-- gdx->rdx for corr index 
    rdx=[pkx[j] for j in gdx]; rdy=[pky[j] for j in gdy]
    tx=npy.take(tx,(gdx),1);ty=npy.take(ty,(gdy),1)
    #print "Applied SVD cut in",round(time.time()-t0,1),'s'

    #--- svd clean
    if int(opt.SVALS)<nbx and int(opt.SVALS)<nby:
        t0=time.time();
        tx,ty=foreach(clean,[tx,ty],threads=2,return_=True)
        print "Cleaned using SVD in",round(time.time()-t0,1),'s'
    else: print "All singulars values retained, no svd clean applied"

    #--- bad bpms to file
    f=open(opt.FILE+'.bad','w')
    print >> f, "@  FILE %s ",opt.FILE
    print >> f, "*  NAME    S    PLANE"
    print >> f, "$   %s     %le   %s "
    for j in range(len(a.H)):
        if j not in rdx:
            print >> f, a.H[j].name, a.H[j].location, "H"
    for j in range(len(a.V)):
        if j not in rdy:
            print >> f, a.V[j].name, a.V[j].location, "V"
    f.close()
    
    #--- good data to file #t0=time.time()
    f = open(opt.FILE+'.new','w')
    f.write('# '+opt.FILE+'\n')
    for j in range(len(gdx)):
        f.write('0 '+a.H[rdx[j]].name+' '+str(a.H[rdx[j]].location)+' ')
        #f.write('%s\n' % ' '.join(['%5.5f' % val for val in \
        #                           a.H[rdx[j]].data]))
        f.write('%s\n' % ' '.join(['%5.5f' % val for val in tx[:,j]]))
    for j in range(len(gdy)):
        f.write('1 '+a.V[rdy[j]].name+' '+str(a.V[rdy[j]].location)+' ')
        #f.write('%s\n' % ' '.join(['%5.5f' % val for val in \
        #                           a.V[rdy[j]].data]))
        f.write('%s\n' % ' '.join(['%5.5f' % val for val in ty[:,j]]))
    f.close();#print "File written in",round(time.time()-t0,1),'s'
    print "Total",round(time.time()-t0,1),'s'


        
if __name__ == "__main__":
    opt,args=parsse();
    

    svdClean()
    
    

