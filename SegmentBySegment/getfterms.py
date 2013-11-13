import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")


from metaclass import twiss
from optparse import OptionParser
import sys


parser = OptionParser()
parser.add_option("-l", "--label",
                help="Label: IP1 IP2 IP3 IP4...",
                metavar="LABEL", default="",dest="label")

parser.add_option("-p", "--path",
                help="Path to model files and output",
                metavar="PATH", default="./", dest="path")

parser.add_option("-f", "--fast",
                help="1 for very fast calculation only writing initvals file",
                metavar="FAST", default="0", dest="fast")

parser.add_option("-e", "--exp",
                help="path to experimental files, only used if fast!=1",
                metavar="EXP", default="./", dest="exp")
parser.add_option("-s", "--start",
                help="Start BPM",
                metavar="START", default="./", dest="start")
parser.add_option("-m", "--method",
                help="Method",
                metavar="ME", default="_free", dest="ME")
parser.add_option("-w", "--w", # Path to Chromaticity functions
                    help="Path to  chromaticity functions, by default this is skiped",
                    metavar="wpath", default="0", dest="wpath")




(options, args) = parser.parse_args()



Model=twiss(options.path+"/twiss_"+options.label+".dat")
method=options.ME
if method=="driven":
	method=""




Model.Cmatrix()

names=Model.NAME

#
# Write initvals.dat with f terms
#
f=open(options.path+"/initvals.dat","w")
f1001=Model.f1001[Model.indx[names[0]]]
f1010=Model.f1010[Model.indx[names[0]]]
print >>f, "f1001r=",f1001.real,";"
print >>f, "f1001i=",f1001.imag,";"
print >>f, "f1010r=",f1010.real,";"
print >>f, "f1010i=",f1010.imag,";"
f.close()




#
# Continue only  if fast option is not 1
# to write all f terms,phase, Disp and W in files
#



if options.fast != "1":


    data=open(options.path+'/sbscouple_'+options.label+'.out','w')
    coupexp=twiss(options.exp+'/getcouple'+method+'.out')
    print options.exp
    ModelPlay=twiss(options.path+"/twiss_"+options.label+"_cor.dat")
    ModelPlay.Cmatrix()
    print >> data,'* NAME   S   f1001 f1001re  f1001im    f1010   f1010re   f1010im  f1001_EXP f1001err_EXP  f1001re_EXP  f1001im_EXP    f1010_EXP  f1010err_EXP   f1010re_EXP   f1010im_EXP     f1001_PLAY   f1001re_PLAY  f1001im_PLAY    f1010_PLAY   f1010re_PLAY   f1010im_PLAY  S_MODEL'
    print >> data,'$  %s    %le   %le    %le  %le   %le    %le    %le  %le   %le    %le  %le   %le    %le    %le   %le  %le   %le    %le    %le   %le  %le %le'

    for name in names:

        f1001=Model.f1001[Model.indx[name]]
        f1010=Model.f1010[Model.indx[name]]
        f1001p=ModelPlay.f1001[ModelPlay.indx[name]]
        f1010p=ModelPlay.f1010[ModelPlay.indx[name]]

        try:
            expi=coupexp.indx[name]
            l=1
        except:
            print "From getfterms.py:",name, "not in experiment"
            l=0
        if l==1:
            print >> data, name, coupexp.S[expi], abs(f1001), f1001.real, f1001.imag, abs(f1010), f1010.real, f1010.imag,\
                  coupexp.F1001W[expi], coupexp.FWSTD1[expi],coupexp.F1001R[expi],coupexp.F1001I[expi],\
                  coupexp.F1010W[expi], coupexp.FWSTD2[expi],coupexp.F1010R[expi],coupexp.F1010I[expi],\
                  abs(f1001p), f1001p.real, f1001p.imag, abs(f1010p), f1010p.real, f1010p.imag,\
                  Model.S[Model.indx[name]]


        #data.write(name+"  "+str(Model.S[Model.indx[name]])+" "+str(abs(f1001))+" "+str(f1001.real)+" "+str(f1001.imag)+" "+str(abs(f1010))+" "+str(f1010.real)+" "+str(f1010.imag)+" \n")


    data.close()



    #
    #  PHASE calculations
    #




    IP=options.label


    t=Model
    tp=ModelPlay
    px=twiss(options.exp+'/getphasetotx'+method+'.out')
    py=twiss(options.exp+'/getphasetoty'+method+'.out')





    def modelIntersect(exp, model):
        bpmsin=[]
        for bpm in exp.NAME:
                try:
                        check=model.indx[bpm.upper()]
                        bpmsin.append([model.S[check], bpm])
                except:
                        print bpm, "Not in Model"
        if len(bpmsin)==0:
                print "Zero intersection of Exp and Model"
                print "Please, provide a good Dictionary"
                print "Now we better leave!"
                sys.exit()
        bpmsin.sort()
        return bpmsin




    bpmsx=modelIntersect(t, px)

    bpmsy=modelIntersect(t, py)


    #    bpmx1=bpmsx[0][1]
    #    bpmy1=bpmsy[0][1]
    bpmx1=options.start
    bpmy1=options.start
    if bpmx1 not in zip(*bpmsx)[1]:   # zip(*a) is like transpose of a list
        print "Selected Start BPM in not in the measurement!"
        print "Quiting getfterms.py"
        sys.exit()

    fx= open(options.path+'/sbsphasext_'+IP+'.out','w')
    print >> fx,"* NAME  S PHASEX  PHASEXT    ERRORX PHASE_PLAY MODEL_S "
    print >> fx,"$ %s    %le  %le  %le        %le     %le       %le  "


    for el in bpmsx:
        tindx=t.indx[el[1]]
        pxindx=px.indx[el[1]]
        mdl=(t.MUX[tindx]-t.MUX[t.indx[bpmx1]] ) % 1
        exp=(px.PHASEX[pxindx]-px.PHASEX[px.indx[bpmx1]] )%1
        exp_mdl=(exp-mdl ) %1
        if exp_mdl >  0.5: exp_mdl=  exp_mdl -1
        mdl_play=-t.MUX[tindx]+tp.MUX[tindx]
        print >>fx,el[1], px.S[pxindx], exp, exp_mdl, px.STDPHX[pxindx], mdl_play, t.S[tindx]

    fx.close()


    fx= open(options.path+'/sbsphaseyt_'+IP+'.out','w')
    print >> fx,"* NAME  S PHASEY  PHASEYT    ERRORY PHASE_PLAY MODEL_S "
    print >> fx,"$ %s   %le  %le   %le        %le    %le     %le "

    for el in bpmsy:
        tindx=t.indx[el[1]]
        pxindx=py.indx[el[1]]
        mdl=(t.MUY[tindx]-t.MUY[t.indx[bpmy1]] ) %1
        exp=(py.PHASEY[pxindx]-py.PHASEY[py.indx[bpmy1]] )%1
        exp_mdl=exp-mdl % 1
        if exp_mdl >  0.5: exp_mdl=  exp_mdl -1
        #if exp_mdl < -0.5: exp_mdl=  exp_mdl +1
        mdl_play=-t.MUY[tindx]+tp.MUY[tindx]
        print >>fx,el[1], py.S[pxindx], exp, exp_mdl, py.STDPHY[pxindx], mdl_play, t.S[tindx]


    fx.close()

#
#   W functions
#

    if options.wpath!="0":
        try:
            wx=twiss(options.wpath+'/wx.out')
            wy=twiss(options.wpath+'/wy.out')
            print "W files loaded correctly"
        except:
            print "No or corrupted W files, skipping"
            options.wpath="0"

    if options.wpath!="0":
        # First WX
        fx= open(options.path+'/sbsWx_'+IP+'.out','w')
        print >> fx,"* NAME  S  WX   WXERR WX_MDL WX_PLAY  PHIX PHIXERR PHIX_MDL PHIX_PLA    MODEL_S "
        print >> fx,"$ %s   %le  %le   %le  %le    %le     %le %le  %le    %le     %le "
        bpmsx=modelIntersect(t, wx)
        for el in bpmsx:
            tindx=t.indx[el[1]]
            tpindx=tp.indx[el[1]]
            wxindx=wx.indx[el[1]]
            print >>fx, el[1], wx.S[wxindx], wx.WX[wxindx], wx.WXERR[wxindx], t.WX[tindx], tp.WX[tpindx], wx.PHIX[wxindx],wx.PHIXERR[wxindx], t.PHIX[tindx], tp.PHIX[tpindx],   t.S[tindx]
        fx.close()
        # Now WY
        fx= open(options.path+'/sbsWy_'+IP+'.out','w')
        print >> fx,"* NAME  S  WY   WYERR WY_MDL WY_PLAY  PHIY PHIYERR PHIY_MDL PHIY_PLAY      MODEL_S "
        print >> fx,"$ %s   %le  %le   %le  %le    %le     %le  %le %le %le %le "
        bpmsx=modelIntersect(t, wy)
        for el in bpmsx:
            tindx=t.indx[el[1]]
            tpindx=tp.indx[el[1]]
            wyindx=wy.indx[el[1]]
            print >>fx, el[1], wy.S[wyindx], wy.WY[wyindx], wy.WYERR[wyindx], t.WY[tindx], tp.WY[tpindx], wy.PHIY[wyindx],wy.PHIYERR[wyindx], t.PHIY[tindx], tp.PHIY[tpindx] ,  t.S[tindx]
        fx.close()




#
# Dispersion
#

    try:
      dx=twiss(options.exp+'/getDx.out')
      ndx=twiss(options.exp+'/getNDx.out')
      dy=twiss(options.exp+'/getDy.out')
      print "All dispersion files loaded"

    except:
      print "No or bad dispersion files in exp dir"
      sys.exit()
    # Writting Dx:
    fx= open(options.path+'/sbsDx_'+IP+'.out','w')
    print >> fx,"* NAME  S  DX   STDDX DX_MDL DX_PLAY MODEL_S "
    print >> fx,"$ %s   %le  %le   %le  %le    %le     %le "
    bpmsx=modelIntersect(t, dx)
    for el in bpmsx:
      tindx=t.indx[el[1]]
      tpindx=tp.indx[el[1]]
      dxindx=dx.indx[el[1]]
      print >>fx, el[1], dx.S[dxindx], dx.DX[dxindx], dx.STDDX[dxindx], t.DX[tindx], tp.DX[tpindx], t.S[tindx]
    fx.close()

    # Writing Dy
    fx= open(options.path+'/sbsDy_'+IP+'.out','w')
    print >> fx,"* NAME  S  DY   STDDY DY_MDL DY_PLAY MODEL_S "
    print >> fx,"$ %s   %le  %le   %le  %le    %le     %le "
    bpmsx=modelIntersect(t, dy)
    for el in bpmsx:
        tindx=t.indx[el[1]]
        tpindx=tp.indx[el[1]]
        dyindx=dy.indx[el[1]]
        print >>fx, el[1], dy.S[dyindx], dy.DY[dyindx], dy.STDDY[dyindx], t.DY[tindx], tp.DY[tpindx], t.S[tindx]
    fx.close()

    # Writing NDx
    fx= open(options.path+'/sbsNDx_'+IP+'.out','w')
    print >> fx,"* NAME  S  NDX   STDNDX NDX_MDL NDX_PLAY MODEL_S "
    print >> fx,"$ %s   %le  %le   %le   %le    %le     %le "
    bpmsx=modelIntersect(t, ndx)
    for el in bpmsx:
      tindx=t.indx[el[1]]
      tpindx=tp.indx[el[1]]
      dxindx=ndx.indx[el[1]]
      print >>fx, el[1], ndx.S[dxindx], ndx.NDX[dxindx], ndx.STDNDX[dxindx], t.DX[tindx]/sqrt(t.BETX[tindx]), tp.DX[tpindx]/sqrt(tp.BETX[tpindx]), t.S[tindx]
    fx.close()











