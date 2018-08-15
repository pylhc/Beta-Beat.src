import os
from math import sqrt, tan, sin, cos, pi

from SegmentBySegment.sbs_writers import sbs_beta_writer
from SegmentBySegment.sbs_writers import sbs_phase_writer
from SegmentBySegment.sbs_writers.sbs_beta_writer import weighted_average_for_SbS_elements
from tfs_files import tfs_file_writer


def _get_ip_tfs_files(save_path):
    file_ip_parabola = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "IP_para.out"))

    file_ip_parabola.add_column_names(["NAME", "S", "BETX", "EBETX", "X", "EX", "BETY", "EBETY", "Y", "EY", "BETX_ph", "EBETX_ph", "BETY_ph", "EBETY_ph"])
    file_ip_parabola.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    file_ip_propagation_x = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "IP_pro_x.out"))

    file_ip_propagation_x.add_column_names(["NAME", "S", "BETSTARX", "EBETSTARX", "X[cm]", "EX[cm]", "BETIPX", "EBETIPX", "ALFIPX", "EALFIPX"])
    file_ip_propagation_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    file_ip_propagation_y = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "IP_pro_y.out"))

    file_ip_propagation_y.add_column_names(["NAME", "S", "BETY", "EBETY", "Y[cm]", "EY[cm]", "BETIPY", "EBETIPY", "ALFIPY", "EALFIPY"])
    file_ip_propagation_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    return file_ip_parabola, file_ip_propagation_x, file_ip_propagation_y


def write_ip(measured_hor_beta_amp, measured_ver_beta_amp,
             beta_x_ip, err_beta_x_ip, alfa_x_ip, err_alfa_x_ip,
             beta_y_ip, err_beta_y_ip, alfa_y_ip, err_alfa_y_ip,
             input_model, measured_hor_phase, measured_ver_phase, element_name,
             accel_inst, save_path):
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
        'input_model': twiss
            twiss input_model
        'measured_hor_phase':
            measured phase advances
        'measured_ver_phase':
            measured phase advances
        'element_name': string
            name of the IP
        'accel': string
            name of the accelerator
        'save_path': string
            where to save file
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''

    file_ip_parabola, file_ip_propagation_x, file_ip_propagation_y = _get_ip_tfs_files(save_path)

    ## constructor
    accel_beam = "B" + str(accel_inst.get_beam())

    ## beta from parabola calculation
    # Map bpm's
    ip_to_bpm_map = {}
    ip_to_bpm_map['IP1'] = ["BPMSW.1L1.", "BPMSW.1R1."]
    ip_to_bpm_map['IP2'] = ["BPMSW.1L2.", "BPMSW.1R2."]
    ip_to_bpm_map['IP5'] = ["BPMSW.1L5.", "BPMSW.1R5."]
    ip_to_bpm_map['IP8'] = ["BPMSW.1L8.", "BPMSW.1R8."]

    if element_name in ip_to_bpm_map:
        # Removed unusedphix, ephix, phiy, ephiy.
        # If needed --> git commit eca47b2d06ab7eab911c830f30dc222a8bcd91a0
        # --vimaier
        (betastar_h, errbx, location_h, errsx, betastar_v, errby, location_v, errsy, betastar_h_p,
         ebetastar_h_p, betastar_v_p, ebetastar_v_p) = _get_ip_from_parabola(ip_to_bpm_map[element_name][0] + accel_beam,
                                                                             ip_to_bpm_map[element_name][1] + accel_beam,
                                                                             measured_hor_beta_amp, measured_ver_beta_amp,
                                                                             measured_hor_phase, measured_ver_phase)
    else:
        betastar_h_p = ebetastar_h_p = betastar_v_p = ebetastar_v_p = betastar_h = errbx = location_h = errsx = betastar_v = errby = location_v = errsy = 0

    file_ip_parabola.add_table_row([element_name, input_model.S[input_model.indx[element_name]],
                                    betastar_h, errbx, location_h, errsx,
                                    betastar_v, errby, location_v, errsy,
                                    betastar_h_p, ebetastar_h_p, betastar_v_p, ebetastar_v_p])
    file_ip_parabola.write_to_file()

    #betstar and waist
    #x
    betastar, errbetastar, waist, errwaist = _get_ip_from_prop(beta_x_ip, err_beta_x_ip, alfa_x_ip, err_alfa_x_ip)
    file_ip_propagation_x.add_table_row([element_name, input_model.S[input_model.indx[element_name]],
                                         betastar, errbetastar, waist, errwaist, round(beta_x_ip, 3), round(err_beta_x_ip, 3), round(alfa_x_ip, 4), round(err_alfa_x_ip, 4)])

    #y
    betastar, errbetastar, waist, errwaist = _get_ip_from_prop(beta_y_ip, err_beta_y_ip, alfa_y_ip, err_alfa_y_ip)
    file_ip_propagation_y.add_table_row([element_name, input_model.S[input_model.indx[element_name]],
                                         betastar, errbetastar, waist, errwaist, round(beta_y_ip, 3), round(err_beta_y_ip, 3), round(alfa_y_ip, 4), round(err_alfa_y_ip, 4)])

    file_ip_propagation_x.write_to_file()
    file_ip_propagation_y.write_to_file()
    ##dispersion


def _get_ip_from_prop(betaip, errbetaip, alfaip, ealfaip):

    #values
    betastar = betaip / (1 + alfaip ** 2)
    waist = alfaip * betaip  # (sign flip)

    #errors
    errwaist = ((ealfaip / abs(alfaip)) + (errbetaip / abs(betaip))) * abs(waist)
    errbetastar = sqrt(errbetaip ** 2 + (((-2 * alfaip) / (1 + alfaip ** 2) ** 2) * alfaip) ** 2)

    waist = waist * 100  # transferring to CM!!!!
    errwaist = errwaist * 100  # transferring to CM!!!!

    return round(betastar, 3), round(errbetastar, 3), round(waist, 3), round(errwaist, 3)


def _get_ip_from_parabola(bpmleft, bpmright, betax, betay, phasex, phasey):
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
        hlp=1
    except:
        hlp=0

    #right horizontal
    try:
        brx=betax.BETX[betax.indx[bpmright]]
        hrp=1
    except:
        hrp=0


    #left vertical
    try:
        bly=betay.BETY[betay.indx[bpmleft]]
        vlp=1
    except:
        vlp=0

    #right vertical
    try:
        bry=betay.BETY[betay.indx[bpmright]]
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


def write_transverse_damper(propagated_models, element_name, input_model, save_path, measured_hor_phase, measured_ver_phase, measured_hor_beta, measured_ver_beta, accel):
    '''
    Function for getting the phase advance between the dampers and pick-ups

    :Parameters:
        'propagated_models.propagation': twiss?
            containing propagation from BPM to damper
        'propagated_models.back_propagation': twiss?
            containing back propagation from BPM to damper
        'element_name': string?
            element name
        'input_model': twiss
            twiss input model
        'save_path': string
            where to save file
        'measured_hor_phase':
            measured phase advances
        'measured_ver_phase':
            measured phase advances
        'errors':
            twissmin,twissminb,twissmax,twissmaxb
    :Return: None
        nothing => writing to file in this function (new/appending)
    '''
    # initial
    damper_map = {}
    bpms_map = {}

    file_transverse_damper = _get_transverse_damper_tfs_file(save_path, accel)

    #beam1
    damper_map['ADTKH.D5L4.B1'] = ['BPMWA.B5L4.B1', 'BPMWA.A5L4.B1']
    damper_map['ADTKH.C5L4.B1'] = ['BPMWA.B5L4.B1', 'BPMWA.A5L4.B1']
    damper_map['ADTKH.B5L4.B1'] = ['BPMWA.B5L4.B1', 'BPMWA.A5L4.B1']
    damper_map['ADTKH.A5L4.B1'] = ['BPMWA.B5L4.B1', 'BPMWA.A5L4.B1']

    damper_map['ADTKV.A5R4.B1'] = ['BPMWA.A5R4.B1', 'BPMWA.B5R4.B1']
    damper_map['ADTKV.B5R4.B1'] = ['BPMWA.A5R4.B1', 'BPMWA.B5R4.B1']
    damper_map['ADTKV.C5R4.B1'] = ['BPMWA.A5R4.B1', 'BPMWA.B5R4.B1']
    damper_map['ADTKV.D5R4.B1'] = ['BPMWA.A5R4.B1', 'BPMWA.B5R4.B1']

    bpms_map['LHCB1_1'] = ['BPMC.9L4.B1', 'BPMC.7L4.B1', 'BPMC.8L4.B1']  # H
    bpms_map['LHCB1_2'] = ['BPMCA.7R4.B1', 'BPMC.9R4.B1', 'BPMC.8R4.B1']  # V

    #beam2
    damper_map['ADTKV.D5L4.B2'] = ['BPMWA.B5L4.B2', 'BPMWA.A5L4.B2']
    damper_map['ADTKV.C5L4.B2'] = ['BPMWA.B5L4.B2', 'BPMWA.A5L4.B2']
    damper_map['ADTKV.B5L4.B2'] = ['BPMWA.B5L4.B2', 'BPMWA.A5L4.B2']
    damper_map['ADTKV.A5L4.B2'] = ['BPMWA.B5L4.B2', 'BPMWA.A5L4.B2']

    damper_map['ADTKH.A5R4.B2'] = ['BPMWA.B5R4.B2', 'BPMWA.A5R4.B2']
    damper_map['ADTKH.B5R4.B2'] = ['BPMWA.B5R4.B2', 'BPMWA.A5R4.B2']
    damper_map['ADTKH.C5R4.B2'] = ['BPMWA.B5R4.B2', 'BPMWA.A5R4.B2']
    damper_map['ADTKH.D5R4.B2'] = ['BPMWA.B5R4.B2', 'BPMWA.A5R4.B2']

    bpms_map['LHCB2_1'] = ['BPMCA.7R4.B2', 'BPMC.9R4.B2', 'BPMC.8R4.B2']  # H
    bpms_map['LHCB2_2'] = ['BPMC.9L4.B2', 'BPMC.7L4.B2', 'BPMC.8L4.B2']  # V

    ## main
    if file_transverse_damper._TfsFileWriter__tfs_table.is_empty():  # To check if the file is empty

        for count in [1, 2]:

            ref_bpm, end_bpm, eight_bpm = bpms_map[accel + "_" + str(count)]

            # hor bpm_pair's
            try:
                in1x = measured_hor_phase.indx[ref_bpm]

            except:
                in1x = -1

            if (in1x < len(measured_hor_phase.indx) - 1) and (in1x != -1):  # when BPM is not last
                if measured_hor_phase.NAME[measured_hor_phase.indx[ref_bpm] + 1] == end_bpm:
                    in2x = 1
                    eightinphasex = 0
                elif measured_hor_phase.NAME[measured_hor_phase.indx[ref_bpm] + 2] == end_bpm:
                    in2x = 1
                    eightinphasex = 1
                else:  # not in measurement
                    in2x = -1
                    eightinphasex = 0
            elif (in1x != -1):  # when BPM is last
                if measured_hor_phase.NAME[0] == end_bpm:
                    in2x = 1
                    eightinphasex = 0
                elif measured_hor_phase.NAME[1] == end_bpm:
                    in2x = 1
                    eightinphasex = 1
                else:  # not in measurement
                    in2x = -1
                    eightinphasex = 0
            else:
                in2x = -1
                eightinphasex = 0

            # ver bpm's
            try:
                in1y = measured_ver_phase.indx[ref_bpm]
            except:
                in1y = -1

            if (in1y < len(measured_ver_phase.indx) - 1) and (in1y != -1):  # when BPM is not last
                if measured_ver_phase.NAME[measured_ver_phase.indx[ref_bpm] + 1] == end_bpm:
                    in2y = 1
                    eightinphasey = 0
                elif measured_ver_phase.NAME[measured_ver_phase.indx[ref_bpm] + 2] == end_bpm:
                    in2y = 1
                    eightinphasey = 1
                else:  # not in measurement
                    in2y = -1
                    eightinphasey = 0
            elif in1y != -1:  # when BPM is last
                if measured_ver_phase.NAME[0] == end_bpm:
                    in2y = 1
                    eightinphasey = 0
                elif measured_ver_phase.NAME[1] == end_bpm:
                    in2y = 1
                    eightinphasey = 1
                else:  # not in measurement
                    in2y = -1
                    eightinphasey = 0
            else:
                in2y = -1
                eightinphasey = 0

            ###### H plane
            if in1x != -1 and in2x != -1:
                if eightinphasex == 1:
                    phaseh = "% .5f" % float(measured_hor_phase.PHASEX[measured_hor_phase.indx[ref_bpm]] + measured_hor_phase.PHASEX[measured_hor_phase.indx[eight_bpm]])
                    errphaseh = "% .5f" % float(measured_hor_phase.STDPHX[measured_hor_phase.indx[ref_bpm]] + measured_hor_phase.STDPHX[measured_hor_phase.indx[eight_bpm]])
                    phasemodelh = "% .5f" % float(measured_hor_phase.PHXMDL[measured_hor_phase.indx[ref_bpm]] + measured_hor_phase.PHXMDL[measured_hor_phase.indx[eight_bpm]])
                else:
                    phaseh = "% .5f" % float(measured_hor_phase.PHASEX[measured_hor_phase.indx[ref_bpm]])
                    errphaseh = "% .5f" % float(measured_hor_phase.STDPHX[measured_hor_phase.indx[ref_bpm]])
                    phasemodelh = "% .5f" % float(measured_hor_phase.PHXMDL[measured_hor_phase.indx[ref_bpm]])
            else:
                print "Horizontal plane not found for transverse dampers pick-ups"
                phaseh = 'NIM'
                errphaseh = 'NIM'
                phasemodelh = 'NIM'

            ###### V plane
            print in1y, in2y, eightinphasey
            if (in1y != -1) and (in2y != -1):
                if eightinphasey == 1:
                    phasev = "% .5f" % float(measured_ver_phase.PHASEY[measured_ver_phase.indx[ref_bpm]] + measured_ver_phase.PHASEY[measured_ver_phase.indx[eight_bpm]])
                    errphasev = "% .5f" % float(measured_ver_phase.STDPHY[measured_ver_phase.indx[ref_bpm]] + measured_ver_phase.STDPHY[measured_ver_phase.indx[eight_bpm]])
                    phasemodelv = "% .5f" % float(measured_ver_phase.PHYMDL[measured_ver_phase.indx[ref_bpm]] + measured_ver_phase.PHYMDL[measured_ver_phase.indx[eight_bpm]])
                else:
                    phasev = "% .5f" % float(measured_ver_phase.PHASEY[measured_ver_phase.indx[ref_bpm]])
                    errphasev = "% .5f" % float(measured_ver_phase.STDPHY[measured_ver_phase.indx[ref_bpm]])
                    phasemodelv = "% .5f" % float(measured_ver_phase.PHYMDL[measured_ver_phase.indx[ref_bpm]])

            else:
                print "Vertical plane not found for transverse dampers pick-ups"
                phasev = 'NIM'
                errphasev = 'NIM'
                phasemodelv = 'NIM'

            file_transverse_damper.add_table_row([ref_bpm, end_bpm, input_model.S[input_model.indx[ref_bpm]], input_model.S[input_model.indx[end_bpm]], phaseh, errphaseh, phasemodelh, phasev, errphasev, phasemodelv])

    if element_name in damper_map:
        bpm_pair = damper_map[element_name]

        mux_bpm1_prop = propagated_models.propagation.MUX[propagated_models.propagation.indx[bpm_pair[0]]]
        mux_bpm2_prop = propagated_models.propagation.MUX[propagated_models.propagation.indx[bpm_pair[1]]]
        mux_elem_prop = propagated_models.propagation.MUX[propagated_models.propagation.indx[element_name]]
        muy_bpm1_prop = propagated_models.propagation.MUY[propagated_models.propagation.indx[bpm_pair[0]]]
        muy_bpm2_prop = propagated_models.propagation.MUY[propagated_models.propagation.indx[bpm_pair[1]]]
        muy_elem_prop = propagated_models.propagation.MUY[propagated_models.propagation.indx[element_name]]

        model_mux_bpm1 = input_model.MUX[input_model.indx[bpm_pair[0]]]
        model_mux_bpm2 = input_model.MUX[input_model.indx[bpm_pair[1]]]
        model_mux_elem = input_model.MUX[input_model.indx[element_name]]

        model_muy_bpm1 = input_model.MUY[input_model.indx[bpm_pair[0]]]
        model_muy_bpm2 = input_model.MUY[input_model.indx[bpm_pair[1]]]
        model_muy_elem = input_model.MUY[input_model.indx[element_name]]

        phase_advance_x_bpm1_prop = abs(mux_elem_prop - mux_bpm1_prop)
        phase_advance_x_bpm2_prop = abs(mux_bpm2_prop - mux_elem_prop)
        phase_advance_y_bpm1_prop = abs(muy_elem_prop - muy_bpm1_prop)
        phase_advance_y_bpm2_prop = abs(muy_bpm2_prop - muy_elem_prop)

        mux_bpm1_back = propagated_models.back_propagation.MUX[propagated_models.back_propagation.indx[bpm_pair[0]]]
        mux_bpm2_back = propagated_models.back_propagation.MUX[propagated_models.back_propagation.indx[bpm_pair[1]]]
        mux_elem_back = propagated_models.back_propagation.MUX[propagated_models.back_propagation.indx[element_name]]
        muy_bpm1_back = propagated_models.back_propagation.MUY[propagated_models.back_propagation.indx[bpm_pair[0]]]
        muy_bpm2_back = propagated_models.back_propagation.MUY[propagated_models.back_propagation.indx[bpm_pair[1]]]
        muy_elem_back = propagated_models.back_propagation.MUY[propagated_models.back_propagation.indx[element_name]]

        phase_advance_x_bpm1_back = abs(mux_bpm1_back - mux_elem_back)
        phase_advance_x_bpm2_back = abs(mux_elem_back - mux_bpm2_back)
        phase_advance_y_bpm1_back = abs(muy_bpm1_back - muy_elem_back)
        phase_advance_y_bpm2_back = abs(muy_elem_back - muy_bpm2_back)

        model_phase_advance_x_bpm1 = abs(model_mux_elem - model_mux_bpm1)
        model_phase_advance_x_bpm2 = abs(model_mux_bpm2 - model_mux_elem)
        model_phase_advance_y_bpm1 = abs(model_muy_elem - model_muy_bpm1)
        model_phase_advance_y_bpm2 = abs(model_muy_bpm2 - model_muy_elem)

        gather_betas_list = [[0, bpm_pair[0]], [0, bpm_pair[1]]]  # The function expects this kind of input
        (beta_x_bpm1, err_beta_x_bpm1, alfa_x_bpm1, err_alfa_x_bpm1,
         beta_x_bpm2, err_beta_x_bpm2, alfa_x_bpm2, err_alfa_x_bpm2) = sbs_beta_writer._get_start_end_betas(gather_betas_list, measured_hor_beta, "X")

        (beta_y_bpm1, err_beta_y_bpm1, alfa_y_bpm1, err_alfa_y_bpm1,
         beta_y_bpm2, err_beta_y_bpm2, alfa_y_bpm2, err_alfa_y_bpm2) = sbs_beta_writer._get_start_end_betas(gather_betas_list, measured_ver_beta, "Y")

        err_phase_x_bpm1_prop = sbs_phase_writer._propagate_error_phase(err_beta_x_bpm1, err_alfa_x_bpm1, phase_advance_x_bpm1_prop, beta_x_bpm1, alfa_x_bpm1)
        err_phase_x_bpm1_back = sbs_phase_writer._propagate_error_phase(err_beta_x_bpm1, err_alfa_x_bpm1, phase_advance_x_bpm1_back, beta_x_bpm1, alfa_x_bpm1)
        err_phase_x_bpm2_prop = sbs_phase_writer._propagate_error_phase(err_beta_x_bpm2, err_alfa_x_bpm2, phase_advance_x_bpm2_prop, beta_x_bpm2, alfa_x_bpm2)
        err_phase_x_bpm2_back = sbs_phase_writer._propagate_error_phase(err_beta_x_bpm2, err_alfa_x_bpm2, phase_advance_x_bpm2_back, beta_x_bpm2, alfa_x_bpm2)
        err_phase_y_bpm1_prop = sbs_phase_writer._propagate_error_phase(err_beta_y_bpm1, err_alfa_y_bpm1, phase_advance_y_bpm1_prop, beta_y_bpm1, alfa_y_bpm1)
        err_phase_y_bpm1_back = sbs_phase_writer._propagate_error_phase(err_beta_y_bpm1, err_alfa_y_bpm1, phase_advance_y_bpm1_back, beta_y_bpm1, alfa_y_bpm1)
        err_phase_y_bpm2_prop = sbs_phase_writer._propagate_error_phase(err_beta_y_bpm2, err_alfa_y_bpm2, phase_advance_y_bpm2_prop, beta_y_bpm2, alfa_y_bpm2)
        err_phase_y_bpm2_back = sbs_phase_writer._propagate_error_phase(err_beta_y_bpm2, err_alfa_y_bpm2, phase_advance_y_bpm2_back, beta_y_bpm2, alfa_y_bpm2)

        # TODO: Is this OK?
        average_phase_advance_x_bpm1, final_error_phase_advance_x_bpm1 = weighted_average_for_SbS_elements(phase_advance_x_bpm1_prop,
                                                                                                           err_phase_x_bpm1_prop,
                                                                                                           phase_advance_x_bpm1_back,
                                                                                                           err_phase_x_bpm1_back)
        average_phase_advance_y_bpm1, final_error_phase_advance_y_bpm1 = weighted_average_for_SbS_elements(phase_advance_y_bpm1_prop,
                                                                                                           err_phase_y_bpm1_prop,
                                                                                                           phase_advance_y_bpm1_back,
                                                                                                           err_phase_y_bpm1_back)
        average_phase_advance_x_bpm2, final_error_phase_advance_x_bpm2 = weighted_average_for_SbS_elements(phase_advance_x_bpm2_prop,
                                                                                                           err_phase_x_bpm2_prop,
                                                                                                           phase_advance_x_bpm2_back,
                                                                                                           err_phase_x_bpm2_back)
        average_phase_advance_y_bpm2, final_error_phase_advance_y_bpm2 = weighted_average_for_SbS_elements(phase_advance_y_bpm2_prop,
                                                                                                           err_phase_y_bpm2_prop,
                                                                                                           phase_advance_y_bpm2_back,
                                                                                                           err_phase_y_bpm2_back)
        file_transverse_damper.add_table_row([bpm_pair[0], element_name,
                                              input_model.S[input_model.indx[bpm_pair[0]]],
                                              input_model.S[input_model.indx[element_name]],
                                              "%.5f" % float(average_phase_advance_x_bpm1),
                                              "%.5f" % float(final_error_phase_advance_x_bpm1),
                                              "%.5f" % float(model_phase_advance_x_bpm1),
                                              "%.5f" % float(average_phase_advance_y_bpm1),
                                              "%.5f" % float(final_error_phase_advance_y_bpm1),
                                              "%.5f" % float(model_phase_advance_y_bpm1)])
        file_transverse_damper.add_table_row([element_name, bpm_pair[1],
                                              input_model.S[input_model.indx[element_name]],
                                              input_model.S[input_model.indx[bpm_pair[1]]],
                                              "%.5f" % float(average_phase_advance_x_bpm2),
                                              "%.5f" % float(final_error_phase_advance_x_bpm2),
                                              "%.5f" % float(model_phase_advance_x_bpm2),
                                              "%.5f" % float(average_phase_advance_y_bpm2),
                                              "%.5f" % float(final_error_phase_advance_y_bpm2),
                                              "%.5f" % float(model_phase_advance_y_bpm2)])
    else:
        print "WARN: Element ", element_name, " Not found in damper_map"

    file_transverse_damper.write_to_file()


def _get_transverse_damper_tfs_file(save_path, accel):
    file_transverse_damper = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "transverse_dampers_phases.out"))

    file_transverse_damper.add_string_descriptor("BEAM", accel)
    file_transverse_damper.add_string_descriptor("NIM", "Not In Measurement")

    file_transverse_damper.add_column_names(["NAME1", "NAME2", "S1", "S2", "PHX", "PHXe", "PHXM", "PHY", "PHYe", "PHYM"])
    file_transverse_damper.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    return file_transverse_damper
