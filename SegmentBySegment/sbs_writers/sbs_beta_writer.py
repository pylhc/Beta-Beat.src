import math
import os

import numpy as np

from tfs_files import tfs_file_writer

import logging

LOGGER = logging.getLogger(__name__)



def write_beta(element_name, is_element, 
               measured_hor_beta, measured_ver_beta,
               input_data, 
               input_model, propagated_models, 
               save_path, beta_summary_file,
               betakind):
    '''
    for beta from phase     and alpha from phase (default): betakind = ""  
    for beta from amplitude and alpha from amplitude:       betakind = "amp"  
    for beta from amplitude and alpha from phase    :       betakind = "bampaphase"  
    for beta from kmod :                                    betakind = "kmod"  
    '''

    file_alfa_x, file_beta_x, file_alfa_y, file_beta_y = _get_beta_tfs_files(element_name, save_path, is_element,betakind)

    model_propagation = propagated_models.propagation
    model_back_propagation = propagated_models.back_propagation
    model_cor = propagated_models.corrected
    model_back_cor = propagated_models.corrected_back_propagation

    if not is_element:
        bpms_list_x = intersect([model_cor, model_propagation, model_back_propagation, model_back_cor, measured_hor_beta])
        bpms_list_y = intersect([model_cor, model_propagation, model_back_propagation, model_back_cor, measured_ver_beta])
    else:
        bpms_list = intersect([model_cor, model_propagation, model_back_propagation, model_back_cor, input_model])
        bpms_list_x = bpms_list_y = bpms_list

#     if betakind=="amp" and (input_data.alphaampXfwd_failed or input_data.alphaampXbak_failed):
#         LOGGER.debug("Alpha calculation failed for plane X, writing zeros to output file")    
#         summary_data_x = _write_zerobeta_for_plane(file_alfa_x, file_beta_x, 
#                           "X", element_name, bpms_list_x, 
#                           input_data,input_model, 
#                           save_path, is_element, beta_summary_file)
#     else:
    summary_data_x = _write_beta_for_plane(file_alfa_x, file_beta_x, "X",
                                       element_name, bpms_list_x, measured_hor_beta,
                                       input_data,input_model,
                                       model_propagation, model_cor, model_back_propagation, model_back_cor,
                                       save_path, is_element, beta_summary_file)

#     if betakind=="amp" and ( input_data.alphaampYfwd_failed or input_data.alphaampYbak_failed):
#         LOGGER.debug("Alpha calculation failed for plane X, writing zeros to output file")    
#         summary_data_y = _write_zerobeta_for_plane(file_alfa_y, file_beta_y, 
#                           "Y", element_name, bpms_list_y, 
#                           input_data,input_model, 
#                           save_path, is_element, beta_summary_file)
#     
#     else:
    summary_data_y = _write_beta_for_plane(file_alfa_y, file_beta_y, "Y",
                                       element_name, bpms_list_y, measured_ver_beta,
                                       input_data,input_model,
                                       model_propagation, model_cor, model_back_propagation, model_back_cor,
                                       save_path, is_element, beta_summary_file)
    
    if is_element:
        _write_summary_data(beta_summary_file, summary_data_x, summary_data_y)
        LOGGER.info("write_beta: is element: returning %f / %f ",summary_data_x[2],summary_data_y[2])
        return (summary_data_x[2], summary_data_x[3], summary_data_x[4], summary_data_x[5],
                summary_data_y[2], summary_data_y[3], summary_data_y[4], summary_data_y[5])
    else:
        LOGGER.info("write_beta: is not element: returning zeros ")
        return [0, 0, 0, 0, 0, 0, 0, 0]


def get_beta_summary_file(save_path):
        beta_summary_file = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbs_summary_bet.out"))
        beta_summary_file.add_column_names(["NAME", "S",
                                            "BETPROPX", "ERRBETPROPX", "ALFPROPX", "ERRALFPROPX",
                                            "BETPROPY", "ERRBETPROPY", "ALFPROPY", "ERRALFPROPY",
                                            "BETXMDL", "BETYMDL", "ALFXMDL", "ALFYMDL", "MDL_S"])
        beta_summary_file.add_column_datatypes(["%s", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le", "%le"])
        return beta_summary_file

def get_betaamp_summary_file(save_path):
        beta_summary_file = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbs_summary_betamp.out"))
        beta_summary_file.add_column_names(["NAME", "S",
                                            "BETPROPX", "ERRBETPROPX", "ALFPROPX", "ERRALFPROPX",
                                            "BETPROPY", "ERRBETPROPY", "ALFPROPY", "ERRALFPROPY",
                                            "BETXMDL", "BETYMDL", "ALFXMDL", "ALFYMDL", "MDL_S"])
        beta_summary_file.add_column_datatypes(["%s", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le", "%le"])
        return beta_summary_file

def get_betaamp_alphaphase_summary_file(save_path):
        beta_summary_file = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbs_summary_betampaphase.out"))
        beta_summary_file.add_column_names(["NAME", "S",
                                            "BETPROPX", "ERRBETPROPX", "ALFPROPX", "ERRALFPROPX",
                                            "BETPROPY", "ERRBETPROPY", "ALFPROPY", "ERRALFPROPY",
                                            "BETXMDL", "BETYMDL", "ALFXMDL", "ALFYMDL", "MDL_S"])
        beta_summary_file.add_column_datatypes(["%s", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le", "%le"])
        return beta_summary_file

def _write_summary_data(beta_summary_file, summary_data_x, summary_data_y):
    
    LOGGER.debug("Length of summary arrays %d %d (need to be >= 9)",len(summary_data_x),len(summary_data_y))
    
    beta_summary_file.add_table_row([summary_data_x[0], summary_data_x[1],
                                     summary_data_x[2], summary_data_x[3], summary_data_x[4], summary_data_x[5],
                                     summary_data_y[2], summary_data_y[3], summary_data_y[4], summary_data_y[5],
                                     summary_data_x[6], summary_data_y[6], summary_data_x[7], summary_data_y[7], summary_data_x[8]])


def _get_beta_tfs_files(element_name, save_path, is_element, betakind):
    '''
    for beta from phase     and alpha from phase (default): betakind = ""  
    for beta from amplitude and alpha from amplitude:       betakind = "amp"  
    for beta from amplitude and alpha from phase    :       betakind = "bampaphase"  
    for beta from kmod :                                    betakind = "kmod"  
    '''
    if betakind is None:
        betakind = "" 
    
    if betakind:
        betakind += "_"
        
    file_beta_x = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetax_" + betakind + element_name + ".out"))
    file_alfa_x = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsalfax_" + betakind + element_name + ".out"))
    file_beta_y = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetay_" + betakind + element_name + ".out"))
    file_alfa_y = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsalfay_" + betakind + element_name + ".out"))

    if not is_element:
        file_beta_x.add_column_names(["NAME", "S", "BETPROPX", "ERRBETPROPX", "BETCORX", "ERRBETCORX", "BETBACKX", "ERRBETBACKX", "BETBACKCORX", "ERRBETBACKCORX", "BETXMDL", "MODEL_S"])
        file_beta_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_x.add_column_names(["NAME", "S", "ALFPROPX", "ERRALFPROPX", "ALFCORX", "ERRALFCORX", "ALFBACKX", "ERRALFBACKX", "ALFBACKCORX", "ERRALFBACKCORX", "ALFXMDL", "MODEL_S"])
        file_alfa_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

        file_beta_y.add_column_names(["NAME", "S", "BETPROPY", "ERRBETPROPY", "BETCORY", "ERRBETCORY", "BETBACKY", "ERRBETBACKY", "BETBACKCORY", "ERRBETBACKCORY", "BETYMDL", "MODEL_S"])
        file_beta_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_y.add_column_names(["NAME", "S", "ALFPROPY", "ERRALFPROPY", "ALFCORY", "ERRALFCORY", "ALFBACKY", "ERRALFBACKY", "ALFBACKCORY", "ERRALFBACKCORY", "ALFYMDL", "MODEL_S"])
        file_alfa_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    else:
        file_beta_x.add_column_names(["NAME", "S", "BETPROPX", "ERRBETPROPX", "BETXMDL", "MODEL_S"])
        file_beta_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_x.add_column_names(["NAME", "S", "ALFPROPX", "ERRALFPROPX", "ALFXMDL", "MODEL_S"])
        file_alfa_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le"])

        file_beta_y.add_column_names(["NAME", "S", "BETPROPY", "ERRBETPROPY", "BETYMDL", "MODEL_S"])
        file_beta_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_y.add_column_names(["NAME", "S", "ALFPROPY", "ERRALFPROPY", "ALFYMDL", "MODEL_S"])
        file_alfa_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le"])

    return file_alfa_x, file_beta_x, file_alfa_y, file_beta_y

def _write_zerobeta_for_plane(file_alfa, file_beta, 
                          plane, element_name, bpms_list, 
                          input_data,input_model, 
                          output_path, is_element, beta_summary_file):
    summary_data = []

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        model_s = input_model.S[input_model.indx[bpm_name]]
        beta_model = getattr(input_model, "BET" + plane)[input_model.indx[bpm_name]]
        alfa_model = getattr(input_model, "ALF" + plane)[input_model.indx[bpm_name]]

        if not is_element:

            file_beta.add_table_row([bpm_name, bpm_s,
                                     0, 0, 0, 0,
                                     0, 0, 0, 0,
                                     beta_model, model_s])
            file_alfa.add_table_row([bpm_name, bpm_s,
                                     0, 0, 0, 0,
                                     0, 0, 0, 0,
                                     alfa_model, model_s])
        else:

            file_alfa.add_table_row([bpm_name, bpm_s, 0, 0, alfa_model, model_s])
            file_beta.add_table_row([bpm_name, bpm_s, 0, 0, beta_model, model_s])
            if element_name in bpm_name:
                summary_data = [bpm_name, bpm_s, 0, 0, 0, 0, beta_model, alfa_model, model_s]

    file_beta.write_to_file()
    file_alfa.write_to_file()
    
    return summary_data


def _write_beta_for_plane(file_alfa, file_beta, 
                          plane, element_name, bpms_list, 
                          measured_beta, 
                          input_data,input_model, 
                          model_propagation, model_cor, model_back_propagation, model_back_cor, 
                          output_path, is_element, beta_summary_file):
    
    LOGGER.debug("_write_beta_for_plane: is_element %r, element_name %r",is_element,element_name)
    
    (beta_start, err_beta_start, alfa_start, err_alfa_start,
     beta_end, err_beta_end, alfa_end, err_alfa_end) = _get_start_end_betas(bpms_list, measured_beta, input_data, plane)
     
    LOGGER.debug("_write_beta_for_plane: %s %s : beta_start %f, err_beta_start %f, alfa_start %f, err_alfa_start %f",
                 element_name, plane,
                 beta_start, err_beta_start, alfa_start, err_alfa_start)
    
    # This can be true only in beta from amplitude
    # it is there the only place where it is set to true
    # so we do not need to ckeck beta kind
    alphaampfwd_failed = getattr(input_data, "alphaamp" + plane + "fwd_failed")
    alphaampbak_failed = getattr(input_data, "alphaamp" + plane + "bak_failed")
        
    
    summary_data = []

    for bpm in bpms_list:
        
        LOGGER.debug("bpm %r",bpm)
        
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        beta_propagation = getattr(model_propagation, "BET" + plane)[model_propagation.indx[bpm_name]]
        beta_back_propagation = getattr(model_back_propagation, "BET" + plane)[model_back_propagation.indx[bpm_name]]

        alfa_propagation = getattr(model_propagation, "ALF" + plane)[model_propagation.indx[bpm_name]]
        alfa_back_propagation = getattr(model_back_propagation, "ALF" + plane)[model_back_propagation.indx[bpm_name]]

        delta_phase_prop = (getattr(model_propagation, "MU" + plane)[model_propagation.indx[bpm_name]]) % 1
        delta_phase_corr = (getattr(model_cor, "MU" + plane)[model_cor.indx[bpm_name]]) % 1
        delta_phase_back = (getattr(model_back_propagation, "MU" + plane)[model_back_propagation.indx[bpm_name]]) % 1
        delta_phase_back_corr = (getattr(model_back_cor, "MU" + plane)[model_back_cor.indx[bpm_name]]) % 1

        err_beta_prop = _propagate_error_beta(err_beta_start, err_alfa_start, delta_phase_prop, beta_propagation, beta_start, alfa_start)
        err_alfa_prop = _propagate_error_alfa(err_beta_start, err_alfa_start, delta_phase_prop, alfa_propagation, beta_start, alfa_start)
        err_beta_back = _propagate_error_beta(err_beta_end, err_alfa_end, delta_phase_back, beta_back_propagation, beta_end, alfa_end)
        err_alfa_back = _propagate_error_beta(err_beta_end, err_alfa_end, delta_phase_back, beta_back_propagation, beta_end, alfa_end)

        if not is_element:
            model_s = measured_beta.S[measured_beta.indx[bpm_name]]
            
            if input_model is None:
                # this the case of matcher, which does not know the model
                # it will now work for bet from amplitude as gatampbetax*.out does not have ALFXMDL columnt 
                beta_model = getattr(measured_beta, "BET" + plane + "MDL")[measured_beta.indx[bpm_name]]
                alfa_model = getattr(measured_beta, "ALF" + plane + "MDL")[measured_beta.indx[bpm_name]]
            else:
                beta_model = getattr(input_model, "BET" + plane )[input_model.indx[bpm_name]]
                alfa_model = getattr(input_model, "ALF" + plane )[input_model.indx[bpm_name]]
            
            
            beta_cor = getattr(model_cor, "BET" + plane)[model_cor.indx[bpm_name]]
            err_beta_cor = _propagate_error_beta(err_beta_start, err_alfa_start, delta_phase_corr, beta_cor, beta_start, alfa_start)

            alfa_cor = getattr(model_cor, "ALF" + plane)[model_cor.indx[bpm_name]]
            err_alfa_cor = _propagate_error_alfa(err_beta_start, err_alfa_start, delta_phase_corr, alfa_cor, beta_start, alfa_start)

            beta_back_cor = getattr(model_back_cor, "BET" + plane)[model_back_cor.indx[bpm_name]]
            err_beta_back_cor = _propagate_error_beta(err_beta_start, err_alfa_start, delta_phase_back_corr, beta_back_cor, beta_start, alfa_start)

            alfa_back_cor = getattr(model_back_cor, "ALF" + plane)[model_back_cor.indx[bpm_name]]
            err_alfa_back_cor = _propagate_error_alfa(err_beta_start, err_alfa_start, delta_phase_back_corr, alfa_back_cor, beta_start, alfa_start)
            
            if alphaampfwd_failed: # it may be true only for beta from amplitude
                beta_propagation = 0
                err_beta_prop    = 0
                beta_cor         = 0
                err_beta_cor     = 0
                alfa_propagation = 0
                err_alfa_prop    = 0
                alfa_cor         = 0
                err_alfa_cor     = 0

            if alphaampbak_failed:
                beta_back_propagation = 0
                err_beta_back         = 0
                beta_back_cor         = 0
                err_beta_back_cor     = 0
                alfa_back_propagation = 0
                err_alfa_back         = 0
                alfa_back_cor         = 0
                err_alfa_back_cor     = 0
                        
            file_beta.add_table_row([bpm_name, bpm_s,
                                     beta_propagation, err_beta_prop, beta_cor, err_beta_cor,
                                     beta_back_propagation, err_beta_back, beta_back_cor, err_beta_back_cor,
                                     beta_model, model_s])
            file_alfa.add_table_row([bpm_name, bpm_s,
                                     alfa_propagation, err_alfa_prop, alfa_cor, err_alfa_cor,
                                     alfa_back_propagation, err_alfa_back, alfa_back_cor, err_alfa_back_cor,
                                     alfa_model, model_s])
        else:
            model_s = input_model.S[input_model.indx[bpm_name]]
            beta_model = getattr(input_model, "BET" + plane)[input_model.indx[bpm_name]]
            alfa_model = getattr(input_model, "ALF" + plane)[input_model.indx[bpm_name]]

            magnet_err_forward = []
            magnet_err_back = []
            if plane == 'X':
                tune = input_model.Q1
            elif plane == 'Y':
                tune = input_model.Q2
            reference_radius = 0.017
            distributed_uncertainty_in_units = 1
            relative_K1L_uncertainty = distributed_uncertainty_in_units * 1E-4 / reference_radius
            
            LOGGER.debug("Using relative_K1L_uncertainty %f ",relative_K1L_uncertainty)
            
            for elem in model_propagation.NAME:
                if abs(getattr(model_propagation, "K1L")[model_propagation.indx[elem]]) > 1E-8: # i.e. if elem == quadrupole
                    K1L_uncertainy = abs(getattr(model_propagation, "K1L")[model_propagation.indx[elem]]) * relative_K1L_uncertainty
                    beta_at_magnet = getattr(model_propagation, "BET" + plane)[model_propagation.indx[elem]]
                    phase_at_magnet = getattr(model_propagation, "MU" + plane)[model_propagation.indx[elem]]
                    phase_at_element = getattr(model_propagation, "MU" + plane)[model_propagation.indx[bpm_name]]
                    delta_phase_from_magnet = phase_at_element - phase_at_magnet
                    if delta_phase_from_magnet > 0:
                        LOGGER.debug("%s : beta_at_magnet %f , delta_phase_from_magnet %f", elem, beta_at_magnet, delta_phase_from_magnet )
                        beta_uncertainty = beta_model * _relative_beta_err_from_magnet_err(K1L_uncertainy, beta_at_magnet, delta_phase_from_magnet, tune)
                        LOGGER.debug("%s : resulting beta_uncertainty %f", elem, beta_uncertainty )
                        magnet_err_forward.append(beta_uncertainty)

            for elem in model_back_propagation.NAME:
                if abs(getattr(model_back_propagation, "K1L")[model_back_propagation.indx[elem]]) > 1E-8: # i.e. if elem == quadrupole
                    K1L_uncertainy = abs(getattr(model_back_propagation, "K1L")[model_back_propagation.indx[elem]]) * relative_K1L_uncertainty
                    beta_at_magnet = getattr(model_back_propagation, "BET" + plane)[model_back_propagation.indx[elem]]
                    phase_at_magnet = getattr(model_back_propagation, "MU" + plane)[model_back_propagation.indx[elem]]
                    phase_at_element = getattr(model_back_propagation, "MU" + plane)[model_back_propagation.indx[bpm_name]]
                    delta_phase_from_magnet = phase_at_element - phase_at_magnet
                    if delta_phase_from_magnet > 0:
                        beta_uncertainty = beta_model * _relative_beta_err_from_magnet_err(K1L_uncertainy, beta_at_magnet, delta_phase_from_magnet, tune)
                        magnet_err_back.append(beta_uncertainty)

            err_beta_prop = np.sqrt(err_beta_prop**2 + sum([x**2 for x in magnet_err_forward])) 
            err_beta_back = np.sqrt(err_beta_back**2 + sum([x**2 for x in magnet_err_back])) 

            
            if alphaampfwd_failed or alphaampbak_failed: # it may be true only for beta from amplitude
                LOGGER.debug("_write_beta_for_plane: %s: beta from amp propagation failed.",element_name)
                averaged_alfa = 0
                final_alfa_error = 0
                averaged_beta = 0
                final_beta_error = 0
            else:
                averaged_beta, final_beta_error = weighted_average_for_SbS_elements(beta_propagation, err_beta_prop, beta_back_propagation, err_beta_back)
                averaged_alfa, final_alfa_error = weighted_average_for_SbS_elements(alfa_propagation, err_alfa_prop, alfa_back_propagation, err_alfa_back)
                

            LOGGER.debug("%s : beta_model %f, averaged_beta %f, beta_propagation %f, beta_back_propagation %f",
                        element_name,beta_model,averaged_beta, beta_propagation,beta_back_propagation)
            
            file_alfa.add_table_row([bpm_name, bpm_s, averaged_alfa, final_alfa_error, alfa_model, model_s])
            file_beta.add_table_row([bpm_name, bpm_s, averaged_beta, final_beta_error, beta_model, model_s])
            if element_name in bpm_name:
                LOGGER.info("_write_beta_for_plane: %s : put to summary averaged_beta %f, final_beta_error %f, beta_model %f ",element_name, averaged_beta, final_beta_error,beta_model)
                summary_data = [bpm_name, bpm_s, averaged_beta, final_beta_error, averaged_alfa, final_alfa_error, beta_model, alfa_model, model_s]
            #else:
            #    LOGGER.info("_write_beta_for_plane: %s : is not in bpm_name list, returning empty summary",element_name)

    file_beta.write_to_file()
    file_alfa.write_to_file()
    return summary_data


def _get_start_end_betas(bpms_list, measured_beta, input_data, plane):
    
    LOGGER.debug("plane %s",plane)
    first_bpm = bpms_list[0][1]

    LOGGER.debug("first_bpm %s",first_bpm)

    beta_start = getattr(measured_beta, "BET" + plane)[measured_beta.indx[first_bpm]]
    if input_data.alphaX_start is None:
        alfa_start = getattr(measured_beta, "ALF" + plane)[measured_beta.indx[first_bpm]]
    else:
        alfa_start = getattr(input_data,'alpha'+plane+'_start')
         

    LOGGER.debug("beta_start = %f ; alfa_start =  %f ;",beta_start, alfa_start)


    #check if BETSTD columns exist (beta amplitude)
    betstd_exists = True
    try:
        getattr(measured_beta, "BET" + plane+"STD")[measured_beta.indx[first_bpm]]
    except:
        betstd_exists = False

    #check if STDBET columns exist
    stdbet_exists = True
    try:
        getattr(measured_beta, "STDBET" + plane)[measured_beta.indx[first_bpm]]
    except:
        stdbet_exists = False
    
    stdalf_exists = True
    try:
        getattr(measured_beta, "STDALF" + plane)[measured_beta.indx[first_bpm]]
    except:
        stdalf_exists = False

    last_bpm = bpms_list[-1][1]

    LOGGER.debug("last_bpm %s",last_bpm)

    beta_end = getattr(measured_beta, "BET" + plane)[measured_beta.indx[last_bpm]]
    
    if input_data.alphaX_end is None:
        alfa_end = -getattr(measured_beta, "ALF" + plane)[measured_beta.indx[last_bpm]]
    else:
        alfa_end = -getattr(input_data,'alpha'+plane+'_end')

    LOGGER.debug("beta_end = %f ; alfa_end =  %f ;",beta_end,alfa_end)
    
    if betstd_exists:
        err_beta_start = getattr(measured_beta, "BET" + plane+"STD")[measured_beta.indx[first_bpm]]
        err_beta_end   = getattr(measured_beta, "BET" + plane+"STD")[measured_beta.indx[last_bpm]]
    elif stdbet_exists:
        err_beta_start = math.sqrt(getattr(measured_beta, "ERRBET" + plane)[measured_beta.indx[first_bpm]] ** 2 + getattr(measured_beta, "STDBET" + plane)[measured_beta.indx[first_bpm]] ** 2)
        err_beta_end   = math.sqrt(getattr(measured_beta, "ERRBET" + plane)[measured_beta.indx[last_bpm]] ** 2 + getattr(measured_beta, "STDBET" + plane)[measured_beta.indx[last_bpm]] ** 2)     
    else:
        err_beta_start = getattr(measured_beta, "ERRBET" + plane)[measured_beta.indx[first_bpm]]
        err_beta_end   = getattr(measured_beta, "ERRBET" + plane)[measured_beta.indx[last_bpm]]
    
    if betstd_exists:
        # beta amplitude does not implement error calculation
        # when it does, it should be passed via input_data
        err_alfa_start = getattr(input_data, "err_alpha" + plane + "_start")
        err_alfa_end   = getattr(input_data, "err_alpha" + plane + "_end")
    elif stdalf_exists:
        err_alfa_start = math.sqrt(getattr(measured_beta, "ERRALF" + plane)[measured_beta.indx[first_bpm]] ** 2 + getattr(measured_beta, "STDALF" + plane)[measured_beta.indx[first_bpm]] ** 2)
        err_alfa_end   = math.sqrt(getattr(measured_beta, "ERRALF" + plane)[measured_beta.indx[last_bpm]] ** 2 + getattr(measured_beta, "STDALF" + plane)[measured_beta.indx[last_bpm]] ** 2)
    else:
        err_alfa_start = getattr(measured_beta, "ERRALF" + plane)[measured_beta.indx[first_bpm]] 
        err_alfa_end   = getattr(measured_beta, "ERRALF" + plane)[measured_beta.indx[last_bpm]]

    return beta_start, err_beta_start, alfa_start, err_alfa_start, beta_end, err_beta_end, alfa_end, err_alfa_end


def _propagate_error_beta(errb0, erra0, dphi, bets, bet0, alf0):
    return math.sqrt((bets*np.sin(4*np.pi*dphi)*alf0/bet0 + bets*np.cos(4*np.pi*dphi)/bet0)**2*errb0**2 + (bets*np.sin(4*np.pi*dphi))**2*erra0**2)  # @IgnorePep8


def _relative_beta_err_from_magnet_err(K1L_uncertainy, beta_at_magnet, delta_phase_from_magnet, tune):
    return K1L_uncertainy * abs(beta_at_magnet * np.sin(4 * np.pi * delta_phase_from_magnet))


def _propagate_error_alfa(errb0, erra0, dphi, alfs, bet0, alf0):
    return math.sqrt(((alfs*((np.sin(4*np.pi*dphi)*alf0/bet0) + (np.cos(4*np.pi*dphi)/bet0))) - (np.cos(4*np.pi*dphi)*alf0/bet0) + (np.sin(4*np.pi*dphi)/bet0))**2*errb0**2 + ((np.cos(4*np.pi*dphi)) - (alfs*np.sin(4*np.pi*dphi)))**2*erra0**2)  # @IgnorePep8


def intersect(list_of_files):
    '''Pure intersection of all bpm names in all files 
     Returns: 
    '''
    if len(list_of_files) == 0:
        raise ValueError("Nothing to intersect!")
    z = list_of_files[0].NAME
    
    for b in list_of_files:
        z = filter(lambda x: x in z and "DRIFT" not in x, b.NAME)
    # SORT by S
    result = []
    x0 = list_of_files[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))
    
    # this makes that the same S is not swapped when sorting 
    result.sort(lambda x, y: 1 if x[0] == y[0] else cmp(x[0], y[0]))
    
    
    if len(result) == 0:
        raise ValueError(
            "The intersection was empty for the files: " +
            str(list_of_files)
        )
    return result


def weighted_average_for_SbS_elements(value1, sigma1, value2, sigma2):
    weighted_average =  (1/sigma1**2 * value1 + 1/sigma2**2 * value2) / (1/sigma1**2 + 1/sigma2**2)  # @IgnorePep8
    uncertainty_of_average = np.sqrt(1 / (1/sigma1**2 + 1/sigma2**2))  # @IgnorePep8
    weighted_rms = np.sqrt(2 * (1/sigma1**2 * (value1 - weighted_average)**2 + 1/sigma2**2 * (value2 - weighted_average)**2) / (1/sigma1**2 + 1/sigma2**2))  # @IgnorePep8
    final_error = np.sqrt(uncertainty_of_average**2 + weighted_rms**2)  # @IgnorePep8
    #print "INPUT OF WeightAV FOR SBS: Value1: " + str(value1)+ " Sigma1: " + str(sigma1)+ " Value 2: " + str(value2) + " Sigma 2: "+str(sigma2)
    return weighted_average, final_error
