


def write_chromatic():
    ##chromatic beta
    if len(chromatic)!=0:

        # horizontal
        chromame=chromatic[0]
        chromamip=chromatic[1]
        chromamap=chromatic[2]

        filex= open(path+'/sbsWx_'+namename+'.out','w')

        if (switch==0):
            bpms=intersect([chromame,chromamip,model,modelp])

            print >> filex,"* NAME  S  WX   WXERR WX_MDL WX_PLAY eWX_play  PHIX PHIXERR PHIX_MDL PHIX_PLAY ePHIX_PLAY    MODEL_S "
            print >> filex,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "
        else:
            bpms=intersect([chromamip,model,modelp])

            print >> filex,"* NAME  S  WX_MDL WX_PLAY eWX_play PHIX_MDL PHIX_PLAY ePHIX_PLAY    MODEL_S "
            print >> filex,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "


        for bpm in bpms:
            s=bpm[0]
            name=bpm[1]

            wmo=model.WX[model.indx[name]]
            phmo=model.PHIX[model.indx[name]]
            smo=model.S[model.indx[name]]

            we=chromamap.WX[chromamap.indx[name]]-chromamip.WX[chromamip.indx[name]]
            phe=chromamap.PHIX[chromamap.indx[name]]-chromamip.PHIX[chromamip.indx[name]]

            if switch==0:
                wme=chromame.WX[chromame.indx[name]]
                ewme=chromame.WXERR[chromame.indx[name]]
                phme=chromame.PHIX[chromame.indx[name]]
                ephme=chromame.PHIXERR[chromame.indx[name]]

                w_cor=modelcor.WX[modelcor.indx[name]]
                ph_cor=modelcor.PHIX[modelcor.indx[name]]

                print >> filex,name,s,wme,ewme,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

            else:

                w_cor=modelp.WX[modelp.indx[name]]
                ph_cor=modelp.PHIX[modelp.indx[name]]

                print >> filex,name,s,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

                if namename in name:

                    filed1=filed1+" "+str(wmo)+" "+str(w_cor)+" "+str(we)+" "+str(phme)+" "+str(ephme)+" "+str(phmo)+" "+str(ph_cor)+" "+str(phe)

        # vertical
        chromame=chromatic[3]
        chromamip=chromatic[4]
        chromamap=chromatic[5]

        filey= open(path+'/sbsWy_'+namename+'.out','w')

        if switch==0:
            bpms=intersect([chromame,chromamip,model,modelp])

            print >> filey,"* NAME  S  WY   WYERR WY_MDL WY_PLAY eWY_play  PHIY PHIYERR PHIYMDL PHIY_PLAY ePHIY_PLAY    MODEL_S "
            print >> filey,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "
        else:
            bpms=intersect([chromamip,model,modelp])

            print >> filey,"* NAME  S  WY_MDL WY_PLAY eWY_play PHIY_MDL PHIY_PLAY ePHIY_PLAY    MODEL_S "
            print >> filey,"$ %s   %le  %le %le   %le  %le  %le    %le     %le %le  %le    %le     %le "

        for bpm in bpms:
            s=bpm[0]
            name=bpm[1]

            wmo=model.WY[model.indx[name]]
            phmo=model.PHIY[model.indx[name]]
            smo=model.S[model.indx[name]]

            we=chromamap.WY[chromamap.indx[name]]-chromamip.WY[chromamip.indx[name]]
            phe=chromamap.PHIY[chromamap.indx[name]]-chromamip.PHIY[chromamip.indx[name]]

            if switch==0:
                wme=chromame.WY[chromame.indx[name]]
                ewme=chromame.WYERR[chromame.indx[name]]
                phme=chromame.PHIY[chromame.indx[name]]
                ephme=chromame.PHIYERR[chromame.indx[name]]

                w_cor=modelcor.WY[modelcor.indx[name]]
                ph_cor=modelcor.PHIY[modelcor.indx[name]]

                print >> filey,name,s,wme,ewme,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

            else:

                w_cor=modelp.WY[modelp.indx[name]]
                ph_cor=modelp.PHIY[modelp.indx[name]]

                print >> filey,name,s,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

                if namename in name:

                    print >>filesum_d,filed1,wmo,w_cor,we,phme,ephme,phmo,ph_cor,phe,smo

        filesum_b.close()
        filesum_c.close()
        filesum_d.close()

        filey.close()