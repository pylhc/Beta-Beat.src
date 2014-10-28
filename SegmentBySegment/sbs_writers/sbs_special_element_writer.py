import os
from math import sqrt, tan, sin, cos, pi


def write_ip(betameA,basetwiss,betatwiss,alfatwiss,model,phasex,phasey,name,accel,path):
    '''
    Function calculating the optics parameters at the IP

    :Parameters:
        'betameA': list
            contaning horizontal and vertical measured amplitude functions
        'basetwiss': list
            list with propogated and back propogated twiss
        'betatwiss': list
            list of twisses for beta errors
        'alfatwiss': list
            list of twisses for alfa errors
        'model': twiss
            twiss model
        'phasex':
            measured phase advances
        'phasey':
            measured phase advances
        'name': string
            name of the IP
        'accel': string
            name of the accelerator
        'save_path': string
            where to save file
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    if os.path.isfile(path+"/IP_para.out"):
        ipfilepara=open(path+"/IP_para.out","a")
    else:
        ipfilepara=open(path+"/IP_para.out","w")
        print >> ipfilepara,"* NAME S BETX EBETX X EX BETY EBETY Y EY BETX_ph EBETX_ph BETY_ph EBETY_ph"
        print >> ipfilepara,"$ %s %le %le %le %le %le %le %le %le %le %le %le %le %le"


    if os.path.isfile(path+"/IP_pro_x.out"):
        ipfileprox=open(path+"/IP_pro_x.out","a")
    else:
        ipfileprox=open(path+"/IP_pro_x.out","w")
        print >> ipfileprox,"* NAME S BETSTARX EBETSTARX X[cm] EX[cm] BETIPX EBETIPX ALFIPX EALFIPX"
        print >> ipfileprox,"$ %s %le %le %le %le %le %le %le %le %le"

    if os.path.isfile(path+"/IP_pro_y.out"):
        ipfileproy=open(path+"/IP_pro_y.out","a")
    else:
        ipfileproy=open(path+"/IP_pro_y.out","w")
        print >> ipfileproy,"* NAME S BETY EBETY Y[cm] EY[cm] BETIPY EBETIPY ALFIPY EALFIPY"
        print >> ipfileproy,"$ %s %le %le %le %le %le %le %le %le %le"

    ## constructor
    if "B1" in accel:
        accelb="B1"
    else:
        accelb="B2"
    twissp=basetwiss[0]
    twissb=basetwiss[1]

    twissbmip=betatwiss[0]
    twissbmib=betatwiss[1]
    twissbmap=betatwiss[2]
    twissbmab=betatwiss[3]

    twissamip=alfatwiss[0]
    twissamib=alfatwiss[1]
    twissamap=alfatwiss[2]
    twissamab=alfatwiss[3]

    ampbetx=betameA[0]
    ampbety=betameA[1]

    ## beta from parabola calculation
    # Map bpm's
    bpmmap={}
    bpmmap['IP1']=["BPMSW.1L1.","BPMSW.1R1."]
    bpmmap['IP2']=["BPMSW.1L2.","BPMSW.1R2."]
    bpmmap['IP5']=["BPMSW.1L5.","BPMSW.1R5."]
    bpmmap['IP8']=["BPMSW.1L8.","BPMSW.1R8."]

    if(name in bpmmap):
        # Removed unusedphix, ephix, phiy, ephiy.
        # If needed --> git commit eca47b2d06ab7eab911c830f30dc222a8bcd91a0
        # --vimaier
        betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_p,ebetastar_h_p,betastar_v_p,ebetastar_v_p=_get_ip_from_para(bpmmap[name][0]+accelb,bpmmap[name][1]+accelb,ampbetx,ampbety,phasex,phasey)
    else:
        betastar_h_p=ebetastar_h_p=betastar_v_p=ebetastar_v_p=betastar_h=errbx=location_h=errsx=betastar_v=errby=location_v=errsy=0

    print >> ipfilepara,name,model.S[model.indx[name]],betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_p,ebetastar_h_p,betastar_v_p,ebetastar_v_p
    ipfilepara.close()

    ## beta from propogations
    #gathering info
    #error alfa

    dalfaxp=1/abs(twissamap.ALFX[twissamap.indx[name]]-twissamip.ALFX[twissamip.indx[name]])
    dalfaxb=1/abs(twissamab.ALFX[twissamab.indx[name]]-twissamib.ALFX[twissamib.indx[name]])

    dalfayp=1/abs(twissamap.ALFY[twissamap.indx[name]]-twissamip.ALFY[twissamip.indx[name]])
    dalfayb=1/abs(twissamab.ALFY[twissamab.indx[name]]-twissamib.ALFY[twissamib.indx[name]])

    #error beta's
    dbxp=1/abs(twissbmap.BETX[twissbmap.indx[name]]-twissbmip.BETX[twissbmip.indx[name]])
    dbxb=1/abs(twissbmab.BETX[twissbmab.indx[name]]-twissbmib.BETX[twissbmib.indx[name]])

    dbyp=1/abs(twissbmap.BETY[twissbmap.indx[name]]-twissbmip.BETY[twissbmip.indx[name]])
    dbyb=1/abs(twissbmab.BETY[twissbmab.indx[name]]-twissbmib.BETY[twissbmib.indx[name]])

    bb=twissb.ALFX[twissb.indx[name]]
    if bb<0.0:
        Balfx=abs(bb)
    else:
        Balfx=-bb

    bb=twissb.ALFY[twissb.indx[name]]
    if bb<0.0:
        Balfy=abs(bb)
    else:
        Balfy=-bb

    alfax=(1/(dalfaxp+dalfaxb))*(dalfaxp*twissp.ALFX[twissp.indx[name]]+Balfx*dalfaxb)

    ealfax=(sqrt((1/dalfaxp)**2+(1/dalfaxb)**2))/2

    alfay=(1/(dalfayp+dalfayb))*(dalfayp*twissp.ALFY[twissp.indx[name]]+Balfy*dalfayb)
    ealfay=(sqrt((1/dalfayp)**2+(1/dalfayb)**2))/2

    betxip=(1/(dbxp+dbxb))*(dbxp*twissp.BETX[twissp.indx[name]]+dbxb*twissb.BETX[twissb.indx[name]])

    betyip=(1/(dbyp+dbyb))*(dbyp*twissp.BETY[twissp.indx[name]]+dbyb*twissb.BETY[twissb.indx[name]])
    errbetxip=sqrt((1/dbxp)**2+(1/dbxb)**2)/2
    errbetyip=sqrt((1/dbyp)**2+(1/dbyb)**2)/2

    #betstar and waist
    #x
    betastar,ebetastar,waist,ewaist=_get_ip_from_prop(betxip,errbetxip,alfax,ealfax)
    print >> ipfileprox,name,model.S[model.indx[name]],betastar,ebetastar,waist,ewaist,round(betxip,3),round(errbetxip,3),round(alfax,4),round(ealfax,4)
    ipfileprox.close()

    #y
    betastar,ebetastar,waist,ewaist=_get_ip_from_prop(betyip,errbetyip,alfay,ealfay)
    print >> ipfileproy,name,model.S[model.indx[name]],betastar,ebetastar,waist,ewaist,round(betyip,3),round(errbetyip,3),round(alfay,4),round(ealfay,4)
    ipfileproy.close()

    ##dispersion


def _get_ip_from_prop(betaip, errbetaip, alfaip, ealfaip):

    #values
    betastar = betaip / (1 + alfaip ** 2)
    waist = alfaip * betaip  # (sign flip)

    #errors
    ewaist=((ealfaip/abs(alfaip))+(errbetaip    /abs(betaip)))*abs(waist)
    ebetastar=sqrt(errbetaip**2+(((-2*alfaip)/(1+alfaip**2)**2)*alfaip)**2 )

    waist=waist*100 # transferring to CM!!!!
    ewaist=ewaist*100 # transferring to CM!!!!

    return round(betastar, 3), round(ebetastar, 3), round(waist, 3), round(ewaist, 3)


def _get_ip_from_para(bpmleft, bpmright, betax, betay, phasex, phasey):
    '''
    Function calculating beta at IP (based on Rioichy thesis)
    b(s)=b*+(s^2/b*)

    :Parameters:
        'bpmleft,bpmright': <dtype?>
            bpm's used to caculate the beta at the IP
        'ampbetx,ampbety': twiss
            twiss containing the horizontal and vertical beta from amp
        'phasex,phasey': <dtype?>
            measured phases
    :Return: <dtype?>
        beta at waist and offset
    '''
    #left horizontal
    try:
        blx=betax.BETX[betax.indx[bpmleft]]
        betax.BETXSTD[betax.indx[bpmleft]]
        hlp=1
    except:
        hlp=0

    #right horizontal
    try:
        brx=betax.BETX[betax.indx[bpmright]]
        betax.BETXSTD[betax.indx[bpmright]]
        hrp=1
    except:
        hrp=0


    #left vertical
    try:
        bly=betay.BETY[betay.indx[bpmleft]]
        betay.BETYSTD[betay.indx[bpmleft]]
        vlp=1
    except:
        vlp=0

    #right vertical
    try:
        bry=betay.BETY[betay.indx[bpmright]]
        betay.BETYSTD[betay.indx[bpmright]]
        vrp=1
    except:
        vrp=0

    #getting phaseadvance
    #inx1=phasex.indx[bpmleft]
    #inx2=phasex.indx[bpmright]
    try:
        inx1=phasex.indx[bpmleft]
        inx2=phasex.indx[bpmright]
        if (inx1==(inx2-1)) or (inx2==0 and inx1==len(phasex.NAME)-1):
            pxp=1
        else:
            pxp=0
    except:
        pxp=0

    try:
        iny1=phasey.indx[bpmleft]
        iny2=phasey.indx[bpmright]
        if iny1==(iny2-1) or (inx2==0 and inx1==len(phasex.NAME)-1):
            pyp=1
        else:
            pyp=0
    except:
        pyp=0


    #horizontal
    if(hrp==1) and (hlp==1) and (pxp==1):

        sxl=phasex.S[inx1]
        sxr=phasex.S[inx2]
        if (sxl==phasex.S[len(phasex.NAME)-1]): #if start point is in center of IP !
            L=sxr
        else:
            L=(sxr-sxl)/2

        phix=phasex.PHASEX[inx1]

        betastar_h=(2*sqrt(blx)*sqrt(brx)*sin(phix*2*pi))/(blx+brx-2*sqrt(blx)*sqrt(brx)*cos(2*pi*phix))*L
        location_h=((blx-brx)/(blx+brx-2*sqrt(blx)*sqrt(brx)*cos(2*pi*phix)))*L



        errbx=0
        errsx=0


        betastar_h_phase=(L)/tan(phix*pi)
        betastar_h_phase_e=(1+cos(2*phix*pi))*L/2


    else:
        print "WARN: Will not calculate horizontal IP"
        phix="NIM"
        betastar_h="NIM"
        location_h="NIM"
        errbx="NIM"
        errsx="NIM"
        betastar_h_phase="NIM"
        betastar_h_phase_e="NIM"
    #vertical
    if (vrp==1) and (vlp==1)and (pyp==1):

        syl=phasey.S[iny1]
        syr=phasey.S[iny2]
        if (syl==phasey.S[len(phasey.NAME)-1]): #if start point is in center of IP !
            L=syr
        else:
            L=(syr-syl)/2


        phiy=phasey.PHASEY[iny1]

        betastar_v=(2*sqrt(bly)*sqrt(bry)*sin(phiy*2*pi))/(bly+bry-2*sqrt(bly)*sqrt(bry)*cos(2*pi*phiy))*L
        location_v=((bly-bry)/(bly+bry-2*sqrt(bly)*sqrt(bry)*cos(2*pi*phiy)))*L

        betastar_v_phase=(L)/tan(phiy*pi)
        betastar_v_phase_e=(1+cos(2*phiy*pi))*L/2

        errby=0
        errsy=0
    else:
        print "WARN: Will not calculate vertical IP"
        phiy="NIM"
        betastar_v="NIM"
        location_v="NIM"
        errby="NIM"
        errsy="NIM"
        betastar_v_phase="NIM"
        betastar_v_phase_e="NIM"

    return betastar_h,errbx,location_h,errsx,betastar_v,errby,location_v,errsy,betastar_h_phase,betastar_h_phase_e,betastar_v_phase,betastar_v_phase_e


def write_transverse_damper(twissp,twissb,element,model,savepath,phasex,phasey,errors):
    '''
    Function for getting the phase advance between the dampers and pick-ups

    :Parameters:
        'twissp': twiss?
            containing propagation from BPM to damper
        'twissb': twiss?
            containing back propagation from BPM to damper
        'element': string?
            element name
        'model': twiss
            twiss model
        'savepath': string
            where to save file
        'phasex':
            measured phase advances
        'phasey':
            measured phase advances
        'errors':
            twissmin,twissminb,twissmax,twissmaxb
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    # initial
    dampermap={}
    bpmmap={}
    path=savepath+"phases.out"
    accel=model.SEQUENCE
    if os.path.exists(path):
        writer=open(path,"a")
        firsttime=0
    else:
        writer=open(path,"w")
        firsttime=1
        print >> writer,"@ BEAM %s",accel
        print >> writer,"@ NIM %s Not In Measurement"
        print >> writer,"* NAME1 NAME2 S1 S2 PHX PHXe PHXM PHY PHYe PHYM"
        print >> writer,"$ %s %s %le %le %le %le %le %le %le %le "


    #beam1
    dampermap['ADTKH.D5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.C5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.B5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']
    dampermap['ADTKH.A5L4.B1']=['BPMWA.B5L4.B1','BPMWA.A5L4.B1']

    dampermap['ADTKV.A5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.B5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.C5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']
    dampermap['ADTKV.D5R4.B1']=['BPMWA.A5R4.B1','BPMWA.B5R4.B1']

    bpmmap['LHCB1_1']=['BPMC.9L4.B1','BPMC.7L4.B1','BPMC.8L4.B1'] #H
    bpmmap['LHCB1_2']=['BPMCA.7R4.B1','BPMC.9R4.B1','BPMC.8R4.B1'] #V

    #beam2
    dampermap['ADTKV.D5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.C5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.B5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']
    dampermap['ADTKV.A5L4.B2']=['BPMWA.B5L4.B2','BPMWA.A5L4.B2']

    dampermap['ADTKH.A5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.B5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.C5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']
    dampermap['ADTKH.D5R4.B2']=['BPMWA.B5R4.B2','BPMWA.A5R4.B2']


    bpmmap['LHCB2_1']=['BPMCA.7R4.B2','BPMC.9R4.B2','BPMC.8R4.B2'] #H
    bpmmap['LHCB2_2']=['BPMC.9L4.B2','BPMC.7L4.B2','BPMC.8L4.B2'] #V

    ## main
    if firsttime==1:

        #
        counter=2
        #

        for count in range(0,counter):

            count=count+1

            print "counter for loop ", count

            bpms=bpmmap[accel+"_"+str(count)]
            refBpm=bpms[0]
            eightbpm=bpms[2]
            endbpm=bpms[1]

            # hor bpms's
            try:
                in1x=phasex.indx[refBpm]

            except:
                in1x=-1

            if (in1x<len(phasex.indx)-1) and (in1x!=-1): # when BPM is not last
                if phasex.NAME[phasex.indx[refBpm]+1]==endbpm:
                    in2x=1
                    eightinphasex=0
                elif phasex.NAME[phasex.indx[refBpm]+2]==endbpm:
                    in2x=1
                    eightinphasex=1
                else: # not in measurement
                    in2x=-1
                    eightinphasex=0
            elif (in1x!=-1): # when BPM is last
                if phasex.NAME[0]==endbpm:
                    in2x=1
                    eightinphasex=0
                elif phasex.NAME[1]==endbpm:
                    in2x=1
                    eightinphasex=1
                else: # not in measurement
                    in2x=-1
                    eightinphasex=0
            else:
                in2x=-1
                eightinphasex=0


            # ver bpm's
            try:
                in1y=phasey.indx[refBpm]
            except:
                in1y=-1

            if (in1y<len(phasey.indx)-1) and (in1y!=-1): # when BPM is not last
                if phasey.NAME[phasey.indx[refBpm]+1]==endbpm:
                    in2y=1
                    eightinphasey=0
                elif phasey.NAME[phasey.indx[refBpm]+2]==endbpm:
                    in2y=1
                    eightinphasey=1
                else: # not in measurement
                    in2y=-1
                    eightinphasey=0
            elif in1y!=-1: # when BPM is last
                if phasey.NAME[0]==endbpm:
                    in2y=1
                    eightinphasey=0
                elif phasey.NAME[1]==endbpm:
                    in2y=1
                    eightinphasey=1
                else: # not in measurement
                    in2y=-1
                    eightinphasey=0
            else:
                in2y=-1
                eightinphasey=0





            ###### H plane
            if (in1x!=-1) and (in2x!=-1):
                if eightinphasex==1:
                    phaseh="%.5f" %float(phasex.PHASEX[phasex.indx[refBpm]]+phasex.PHASEX[phasex.indx[eightbpm]])
                    errphaseh="%.5f" %float(phasex.STDPHX[phasex.indx[refBpm]]+phasex.STDPHX[phasex.indx[eightbpm]])
                    phasemodelh="%.5f" %float(phasex.PHXMDL[phasex.indx[refBpm]]+phasex.PHXMDL[phasex.indx[eightbpm]])
                else:
                    phaseh="%.5f" %float(phasex.PHASEX[phasex.indx[refBpm]])
                    errphaseh="%.5f" %float(phasex.STDPHX[phasex.indx[refBpm]])
                    phasemodelh="%.5f" %float(phasex.PHXMDL[phasex.indx[refBpm]])
            else:
                print "Horizontal plane not found for transverse dampers pick-ups"
                phaseh='NIM'
                errphaseh='NIM'
                phasemodelh='NIM'

            ###### V plane
            print in1y,in2y,eightinphasey
            if (in1y!=-1) and (in2y!=-1):
                if eightinphasey==1:
                    phasev="%.5f" %float(phasey.PHASEY[phasey.indx[refBpm]]+phasey.PHASEY[phasey.indx[eightbpm]])
                    errphasev="%.5f" %float(phasey.STDPHY[phasey.indx[refBpm]]+phasey.STDPHY[phasey.indx[eightbpm]])
                    phasemodelv="%.5f" %float(phasey.PHYMDL[phasey.indx[refBpm]]+phasey.PHYMDL[phasey.indx[eightbpm]])
                else:

                    phasev="%.5f" %float(phasey.PHASEY[phasey.indx[refBpm]])
                    errphasev="%.5f" %float(phasey.STDPHY[phasey.indx[refBpm]])
                    phasemodelv="%.5f" %float(phasey.PHYMDL[phasey.indx[refBpm]])

            else:
                print "Vertical plane not found for transverse dampers pick-ups"
                phasev='NIM'
                errphasev='NIM'
                phasemodelv='NIM'


            print >> writer,refBpm,endbpm,model.S[model.indx[refBpm]],model.S[model.indx[endbpm]],phaseh,errphaseh,phasemodelh,phasev,errphasev,phasemodelv

    if element in dampermap:
        bpms=dampermap[element]
        passs=1
    else:
        print "WARN: Element ",element," Not found in dampermap"
        passs=0

    if passs==1:
        muxbpm1=twissp.MUX[twissp.indx[bpms[0]]]
        muxbpm2=twissp.MUX[twissp.indx[bpms[1]]]
        mux=twissp.MUX[twissp.indx[element]]
        muybpm1=twissp.MUY[twissp.indx[bpms[0]]]
        muybpm2=twissp.MUY[twissp.indx[bpms[1]]]
        muy=twissp.MUY[twissp.indx[element]]

        muxbpm1m=model.MUX[model.indx[bpms[0]]]
        muxbpm2m=model.MUX[model.indx[bpms[1]]]
        muxm=model.MUX[model.indx[element]]

        muybpm1m=model.MUY[model.indx[bpms[0]]]
        muybpm2m=model.MUY[model.indx[bpms[1]]]
        muym=model.MUY[model.indx[element]]

        pha1xp=abs(mux-muxbpm1)
        pha2xp=abs(muxbpm2-mux)
        pha1yp=abs(muy-muybpm1)
        pha2yp=abs(muybpm2-muy)

        muxbpm1=twissb.MUX[twissb.indx[bpms[0]]]
        muxbpm2=twissb.MUX[twissb.indx[bpms[1]]]
        mux=twissb.MUX[twissb.indx[element]]
        muybpm1=twissb.MUY[twissb.indx[bpms[0]]]
        muybpm2=twissb.MUY[twissb.indx[bpms[1]]]
        muy=twissb.MUY[twissb.indx[element]]


        pha1xb=abs(muxbpm1-mux)
        pha2xb=abs(mux-muxbpm2)
        pha1yb=abs(muybpm1-muy)
        pha2yb=abs(muy-muybpm2)


        pha1xm=abs(muxm-muxbpm1m)
        pha2xm=abs(muxbpm2m-muxm)
        pha1ym=abs(muym-muybpm1m)
        pha2ym=abs(muybpm2m-muym)




        #mean
        pha1x=(pha1xp+pha1xb)/2
        pha1y=(pha1yp+pha1yb)/2
        pha2x=(pha2xp+pha2xb)/2
        pha2y=(pha2yp+pha2yb)/2

        #error
        twissmin=errors[0]
        twissminb=errors[1]
        twissmax=errors[2]
        twissmaxb=errors[3]

        minn=abs(abs(twissmin.MUX[twissmin.indx[element]]-twissmin.MUX[twissmin.indx[bpms[0]]])-abs(twissminb.MUX[twissminb.indx[element]]-twissminb.MUX[twissminb.indx[bpms[0]]]))
        maxx=abs(abs(twissmax.MUX[twissmax.indx[element]]-twissmax.MUX[twissmax.indx[bpms[0]]])-abs(twissmaxb.MUX[twissmaxb.indx[element]]-twissmaxb.MUX[twissmaxb.indx[bpms[0]]]))
        pha1xe=minn+maxx


        minn=abs(abs(twissmin.MUX[twissmin.indx[bpms[1]]]-twissmin.MUX[twissmin.indx[element]])-abs(twissminb.MUX[twissminb.indx[bpms[1]]]-twissminb.MUX[twissminb.indx[element]]))
        maxx=abs(abs(twissmax.MUX[twissmax.indx[bpms[1]]]-twissmax.MUX[twissmax.indx[element]])-abs(twissmaxb.MUX[twissmaxb.indx[bpms[1]]]-twissmaxb.MUX[twissmaxb.indx[element]]))
        pha2xe=minn+maxx

        minn=abs(abs(twissmin.MUY[twissmin.indx[element]]-twissmin.MUY[twissmin.indx[bpms[0]]])-abs(twissminb.MUY[twissminb.indx[element]]-twissminb.MUY[twissminb.indx[bpms[0]]]))
        maxx=abs(abs(twissmax.MUY[twissmax.indx[element]]-twissmax.MUY[twissmax.indx[bpms[0]]])-abs(twissmaxb.MUY[twissmaxb.indx[element]]-twissmaxb.MUY[twissmaxb.indx[bpms[0]]]))
        pha1ye=minn+maxx

        minn=abs(abs(twissmin.MUY[twissmin.indx[bpms[1]]]-twissmin.MUY[twissmin.indx[element]])-abs(twissminb.MUY[twissminb.indx[bpms[1]]]-twissminb.MUY[twissminb.indx[element]]))
        maxx=abs(abs(twissmax.MUY[twissmax.indx[bpms[1]]]-twissmax.MUY[twissmax.indx[element]])-abs(twissmaxb.MUY[twissmaxb.indx[bpms[1]]]-twissmaxb.MUY[twissmaxb.indx[element]]))
        pha2ye=minn+maxx

        print >> writer,bpms[0],element,model.S[model.indx[bpms[0]]],model.S[model.indx[element]],"%.5f" %float(pha1x),"%.5f" %float(pha1xe),"%.5f" %float(pha1xm),"%.5f" %float(pha1y),"%.5f" %float(pha1ye),"%.5f" %float(pha1ym)
        print >> writer,element,bpms[1],model.S[model.indx[element]],model.S[model.indx[bpms[1]]],"%.5f" %float(pha2x),"%.5f" %float(pha2xe),"%.5f" %float(pha2xm),"%.5f" %float(pha2y),"%.5f" %float(pha2ye),"%.5f" %float(pha2ym)

    writer.close()
