import os
import numpy as np
import math
from Utilities import tfs_file_writer
import SegmentBySegment


def write_coupling(element_name, is_element, measured_coupling, input_model, propagated_models, save_path):

    file_f_terms, file_C_terms = _get_coupling_tfs_files(element_name, save_path, is_element)

    model_propagation = propagated_models.propagation
    model_back_propagation = propagated_models.back_propagation
    model_cor = propagated_models.corrected
    model_back_cor = propagated_models.corrected_back_propagation

    model_propagation.Cmatrix()
    model_back_propagation.Cmatrix()
    model_cor.Cmatrix()

    bpms_list = SegmentBySegment.intersect([measured_coupling, model_propagation])

    _write_f_terms(file_f_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor)

    _write_C_terms(file_C_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor)


def _get_coupling_tfs_files(element_name, save_path, is_element):
    file_f_terms = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbscouple_" + element_name + ".out"))  # TODO: Probably we want a better name for these files
    file_C_terms = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbscouple2_" + element_name + ".out"))  # TODO: Probably we want a better name for these files

    if not is_element:
        file_f_terms.add_column_names(["NAME", "S",
                                       "F1001ABSPROP", "ERRF1001ABSPROP", "F1001REPROP", "ERRF1001REPROP", "F1001IMPROP", "ERRF1001IMPROP",
                                       "F1010ABSPROP", "ERRF1010ABSPROP", "F1010REPROP", "ERRF1010REPROP", "F1010IMPROP", "ERRF1010IMPROP",
                                       "F1001ABSCOR", "ERRF1001ABSCOR", "F1001RECOR", "ERRF1001RECOR", "F1001IMCOR", "ERRF1001IMCOR",
                                       "F1010ABSCOR", "ERRF1010ABSCOR", "F1010RECOR", "ERRF1010RECOR", "F1010IMCOR", "ERRF1010IMCOR",
                                       "F1001ABSBACK", "ERRF1001ABSBACK", "F1001REBACK", "ERRF1001REBACK", "F1001IMBACK", "ERRF1001IMBACK",
                                       "F1010ABSBACK", "ERRF1010ABSBACK", "F1010REBACK", "ERRF1010REBACK", "F1010IMBACK", "ERRF1010IMBACK",
                                       "F1001ABSBACKCOR", "ERRF1001ABSBACKCOR", "F1001REBACKCOR", "ERRF1001REBACKCOR", "F1001IMBACKCOR", "ERRF1001IMBACKCOR",
                                       "F1010ABSBACKCOR", "ERRF1010ABSBACKCOR", "F1010REBACKCOR", "ERRF1010REBACKCOR", "F1010IMBACKCOR", "ERRF1010IMBACKCOR",
                                       "MODEL_F1001ABS", "MODEL_F1001RE", "MODEL_F1001IM", "MODEL_F1010ABS", "MODEL_F1010RE", "MODEL_F1010IM", "MODEL_S"])
        file_f_terms.add_column_datatypes(["%bpm_s", "%le", 
                                           "%le", "%le", "%le", "%le", "%le", "%le", 
                                           "%le", "%le", "%le", "%le", "%le", "%le", 
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

        file_C_terms.add_column_names(["NAME", "S", "C11Mo", "C12Mo", "C21Mo", "C22Mo", "ANDMo", "C11_cor", "eC11_cor", "C12_cor", "eC12_cor", "C21_cor", "eC21_cor", "C22_cor", "eC22_cor", "ANG_cor", "eANG_cor"])
        file_C_terms.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    else:
        file_f_terms.add_column_names(["NAME", "S",
                                       "F1001ABSPROP", "ERRF1001ABSPROP", "F1001REPROP", "ERRF1001REPROP", "F1001IMPROP", "ERRF1001IMPROP",
                                       "F1010ABSPROP", "ERRF1010ABSPROP", "F1010REPROP", "ERRF1010REPROP", "F1010IMPROP", "ERRF1010IMPROP",
                                       "MODEL_F1001ABS", "MODEL_F1001RE", "MODEL_F1001IM", "MODEL_F1010ABS", "MODEL_F1010RE", "MODEL_F1010IM", "MODEL_S"])
        file_f_terms.add_column_datatypes(["%bpm_s", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    return file_f_terms, file_C_terms


def _write_f_terms(file_f_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor):

    first_bpm = bpms_list[0][1]
    last_bpm = bpms_list[-1][1]

    (f1001ab_ini, f1001_std_ini, p1001_ini, p1001_std_ini,
     f1001ab_end, f1001_std_end, p1001_end, p1001_std_end) = _get_start_end_f1001(measured_coupling, first_bpm, last_bpm)
    (f1010ab_ini, f1010_std_ini, p1010_ini, p1010_std_ini,
     f1010ab_end, f1010_std_end, p1010_end, p1010_std_end) = _get_start_end_f1010(measured_coupling, first_bpm, last_bpm)

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        model_s = input_model.S[input_model.indx[bpm_name]]
        model_f1001 = input_model.f1001[input_model.indx[bpm_name]]
        model_f1010 = input_model.f1010[input_model.indx[bpm_name]]

        delta_phase_prop_x = model_propagation.MUX[model_propagation.indx[bpm_name]] % 1
        delta_phase_prop_y = model_propagation.MUY[model_propagation.indx[bpm_name]] % 1

        f1001_prop = model_propagation.f1001[model_propagation.indx[bpm_name]]
        err_f1001re_prop = _propagate_error_coupling_1001_re(f1001ab_ini, p1001_ini, delta_phase_prop_x, delta_phase_prop_y, f1001_std_ini, p1001_std_ini)
        err_f1001im_prop = _propagate_error_coupling_1001_im(f1001ab_ini, p1001_ini, delta_phase_prop_x, delta_phase_prop_y, f1001_std_ini, p1001_std_ini)
        err_f1001abs_prop = f1001ab_ini
        f1010_prop = model_propagation.f1010[model_propagation.indx[bpm_name]]
        err_f1010re_prop = _propagate_error_coupling_1010_re(f1010ab_ini, p1010_ini, delta_phase_prop_x, delta_phase_prop_y, f1010_std_ini, p1010_std_ini)
        err_f1010im_prop = _propagate_error_coupling_1010_im(f1010ab_ini, p1010_ini, delta_phase_prop_x, delta_phase_prop_y, f1010_std_ini, p1010_std_ini)
        err_f1010abs_prop = f1010ab_ini

        delta_phase_back_x = model_back_propagation.MUX[model_back_propagation.indx[bpm_name]] % 1
        delta_phase_back_y = model_back_propagation.MUY[model_back_propagation.indx[bpm_name]] % 1

        f1001_back = model_back_propagation.f1001[model_back_propagation.indx[bpm_name]]
        err_f1001re_back = _propagate_error_coupling_1001_re(f1001ab_end, p1001_end, delta_phase_back_x, delta_phase_back_y, f1001_std_end, p1001_std_end)
        err_f1001im_back = _propagate_error_coupling_1001_im(f1001ab_end, p1001_end, delta_phase_back_x, delta_phase_back_y, f1001_std_end, p1001_std_end)
        err_f1001abs_back = f1001ab_end
        f1010_back = model_back_propagation.f1010[model_back_propagation.indx[bpm_name]]
        err_f1010re_back = _propagate_error_coupling_1010_re(f1010ab_end, p1010_end, delta_phase_back_x, delta_phase_back_y, f1010_std_end, p1010_std_end)
        err_f1010im_back = _propagate_error_coupling_1010_im(f1010ab_end, p1010_end, delta_phase_back_x, delta_phase_back_y, f1010_std_end, p1010_std_end)
        err_f1010abs_back = f1010ab_end

        if not is_element:
            delta_phase_corr_x = model_cor.MUX[model_cor.indx[bpm_name]] % 1
            delta_phase_corr_y = model_cor.MUY[model_cor.indx[bpm_name]] % 1

            f1001_corr = model_cor.f1001[model_cor.indx[bpm_name]]
            err_f1001re_corr = _propagate_error_coupling_1001_re(f1001ab_ini, p1001_ini, delta_phase_corr_x, delta_phase_corr_y, f1001_std_ini, p1001_std_ini)
            err_f1001im_corr = _propagate_error_coupling_1001_im(f1001ab_ini, p1001_ini, delta_phase_corr_x, delta_phase_corr_y, f1001_std_ini, p1001_std_ini)
            err_f1001abs_corr = f1001ab_ini
            f1010_corr = model_cor.f1010[model_cor.indx[bpm_name]]
            err_f1010re_corr = _propagate_error_coupling_1010_re(f1010ab_ini, p1010_ini, delta_phase_corr_x, delta_phase_corr_y, f1010_std_ini, p1010_std_ini)
            err_f1010im_corr = _propagate_error_coupling_1010_im(f1010ab_ini, p1010_ini, delta_phase_corr_x, delta_phase_corr_y, f1010_std_ini, p1010_std_ini)
            err_f1010abs_corr = f1010ab_ini

            delta_phase_back_corr_x = model_back_cor.MUX[model_back_cor.indx[bpm_name]] % 1
            delta_phase_back_corr_y = model_back_cor.MUY[model_back_cor.indx[bpm_name]] % 1

            f1001_back_corr = model_back_cor.f1001[model_back_cor.indx[bpm_name]]
            err_f1001re_back_corr = _propagate_error_coupling_1001_re(f1001ab_end, p1001_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1001_std_end, p1001_std_end)
            err_f1001im_back_corr = _propagate_error_coupling_1001_im(f1001ab_end, p1001_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1001_std_end, p1001_std_end)
            err_f1001abs_back_corr = f1001ab_end
            f1010_back_corr = model_back_cor.f1010[model_back_cor.indx[bpm_name]]
            err_f1010re_back_corr = _propagate_error_coupling_1010_re(f1010ab_end, p1010_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1010_std_end, p1010_std_end)
            err_f1010im_back_corr = _propagate_error_coupling_1010_im(f1010ab_end, p1010_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1010_std_end, p1010_std_end)
            err_f1010abs_back_corr = f1010ab_end

            file_f_terms.add_table_row([bpm_name, bpm_s,
                                        abs(f1001_prop), err_f1001abs_prop, f1001_prop.real, err_f1001re_prop, f1001_prop.imag, err_f1001im_prop,
                                        abs(f1010_prop), err_f1010abs_prop, f1010_prop.real, err_f1010re_prop, f1010_prop.imag, err_f1010im_prop,
                                        abs(f1001_corr), err_f1001abs_corr, f1001_corr.real, err_f1001re_corr, f1001_corr.imag, err_f1001im_corr,
                                        abs(f1010_corr), err_f1010abs_corr, f1010_corr.real, err_f1010re_corr, f1010_corr.imag, err_f1010im_corr,
                                        abs(f1001_back), err_f1001abs_back, f1001_back.real, err_f1001re_back, f1001_back.imag, err_f1001im_back,
                                        abs(f1010_back), err_f1010abs_back, f1010_back.real, err_f1010re_back, f1010_back.imag, err_f1010im_back,
                                        abs(f1001_back_corr), err_f1001abs_back_corr, f1001_back_corr.real, err_f1001re_back_corr, f1001_back_corr.imag, err_f1001im_back_corr,
                                        abs(f1010_back_corr), err_f1010abs_back_corr, f1010_back_corr.real, err_f1010re_back_corr, f1010_back_corr.imag, err_f1010im_back_corr,
                                        abs(model_f1001), model_f1001.real, model_f1001.imag, abs(model_f1010), model_f1010.real, model_f1010.imag, model_s])
        else:
            average_f1001ab, final_f1001abs_err = SegmentBySegment.weighted_average_for_SbS_elements(abs(f1001_prop), err_f1001abs_prop, abs(f1001_back), err_f1001abs_back)
            average_f1001re, final_f1001re_err = SegmentBySegment.weighted_average_for_SbS_elements(f1001_prop.real, err_f1001re_prop, f1001_back.real, err_f1001re_back)
            average_f1001im, final_f1001im_err = SegmentBySegment.weighted_average_for_SbS_elements(f1001_prop.imag, err_f1001im_prop, f1001_back.imag, err_f1001im_back)

            average_f1010ab, final_f1010abs_err = SegmentBySegment.weighted_average_for_SbS_elements(abs(f1010_prop), err_f1010abs_prop, abs(f1010_back), err_f1010abs_back)
            average_f1010re, final_f1010re_err = SegmentBySegment.weighted_average_for_SbS_elements(f1010_prop.real, err_f1010re_prop, f1010_back.real, err_f1010re_back)
            average_f1010im, final_f1010im_err = SegmentBySegment.weighted_average_for_SbS_elements(f1010_prop.imag, err_f1010im_prop, f1010_back.imag, err_f1010im_back)

            file_f_terms.add_table_row([bpm_name, bpm_s,
                                        average_f1001ab, final_f1001abs_err, average_f1001re, final_f1001re_err, average_f1001im, final_f1001im_err,
                                        average_f1010ab, final_f1010abs_err, average_f1010re, final_f1010re_err, average_f1010im, final_f1010im_err,
                                        abs(model_f1001), model_f1001.real, model_f1001.imag, abs(model_f1010), model_f1010.real, model_f1010.imag, model_s])

    file_f_terms.write_to_file()


def _write_C_terms(file_C_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor):

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        # TODO: All this
        if not is_element:
            file_C_terms.add_table_row([])
        else:
            file_C_terms.add_table_row([])


def _get_start_end_f1001(measured_coupling, first_bpm, last_bpm):
    f1001ab_ini = measured_coupling.F1001W[measured_coupling.indx[first_bpm]]
    f1001_std_ini = measured_coupling.FWSTD1[measured_coupling.indx[first_bpm]]
    p1001_ini = measured_coupling.Q1001[measured_coupling.indx[first_bpm]]
    p1001_std_ini = measured_coupling.Q1001STD[measured_coupling.indx[first_bpm]]

    f1001ab_end = measured_coupling.F1001W[measured_coupling.indx[last_bpm]]
    f1001_std_end = measured_coupling.FWSTD1[measured_coupling.indx[last_bpm]]
    p1001_end = measured_coupling.Q1001[measured_coupling.indx[last_bpm]]
    p1001_std_end = measured_coupling.Q1001STD[measured_coupling.indx[last_bpm]]

    return f1001ab_ini, f1001_std_ini, p1001_ini, p1001_std_ini, f1001ab_end, f1001_std_end, p1001_end, p1001_std_end


def _get_start_end_f1010(measured_coupling, first_bpm, last_bpm):
    f1010ab_ini = measured_coupling.F1010W[measured_coupling.indx[first_bpm]]
    f1010_std_ini = measured_coupling.FWSTD1[measured_coupling.indx[first_bpm]]
    p1010_ini = measured_coupling.Q1010[measured_coupling.indx[first_bpm]]
    p1010_std_ini = measured_coupling.Q1010STD[measured_coupling.indx[first_bpm]]

    f1010ab_end = measured_coupling.F1010W[measured_coupling.indx[last_bpm]]
    f1010_std_end = measured_coupling.FWSTD2[measured_coupling.indx[last_bpm]]
    p1010_end = measured_coupling.Q1010[measured_coupling.indx[last_bpm]]
    p1010_std_end = measured_coupling.Q1010STD[measured_coupling.indx[last_bpm]]

    return f1010ab_ini, f1010_std_ini, p1010_ini, p1010_std_ini, f1010ab_end, f1010_std_end, p1010_end, p1010_std_end


def _propagate_error_coupling_1001_re(f1001ab_ini, p1001_ini, phasex, phasey, f1001_std_ini, p1001_std_ini):
    math.sqrt((f1001_std_ini * np.cos(p1001_ini - phasex + phasey))**2 + (p1001_std_ini * f1001ab_ini * np.sin(p1001_ini - phasex + phasey))**2)  # @IgnorePep8


def _propagate_error_coupling_1001_im(f1001ab_ini, p1001_ini, phasex, phasey, f1001_std_ini, p1001_std_ini):
    math.sqrt((f1001_std_ini * np.sin(p1001_ini - phasex + phasey))**2 + (p1001_std_ini * f1001ab_ini * np.cos(p1001_ini - phasex + phasey))**2)  # @IgnorePep8


def _propagate_error_coupling_1010_re(f1010ab_ini, p1010_ini, phasex, phasey, f1010_std_ini, p1010_std_ini):
    math.sqrt((f1010_std_ini * np.cos(p1010_ini - phasex - phasey))**2 + (p1010_std_ini * f1010ab_ini * np.sin(p1010_ini - phasex - phasey))**2)  # @IgnorePep8


def _propagate_error_coupling_1010_im(f1010ab_ini, p1010_ini, phasex, phasey, f1010_std_ini, p1010_std_ini):
    math.sqrt((f1010_std_ini * np.sin(p1010_ini - phasex - phasey))**2 + (p1010_std_ini * f1010ab_ini * np.cos(p1010_ini - phasex - phasey))**2)  # @IgnorePep8
