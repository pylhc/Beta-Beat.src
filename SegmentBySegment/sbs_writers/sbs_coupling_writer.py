import math
import os

import numpy as np

from sbs_beta_writer import intersect, weighted_average_for_SbS_elements
from tfs_files import tfs_file_writer


def write_coupling(element_name, is_element, measured_coupling, input_model, propagated_models, save_path, coupling_summary_file):

    file_f_terms, file_C_terms = _get_coupling_tfs_files(element_name, save_path, is_element)

    model_propagation = propagated_models.propagation
    model_back_propagation = propagated_models.back_propagation
    model_cor = propagated_models.corrected
    model_back_cor = propagated_models.corrected_back_propagation

    model_propagation.Cmatrix()
    model_back_propagation.Cmatrix()
    model_cor.Cmatrix()
    model_back_cor.Cmatrix()

    if not is_element:
        bpms_list = intersect([model_cor, model_propagation, model_back_propagation, model_back_cor, measured_coupling])
    else:
        bpms_list = intersect([model_cor, model_propagation, model_back_propagation, model_back_cor, input_model])

    summary_data_f = _write_f_terms(file_f_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor)

    summary_data_C = _write_C_terms(file_C_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor)

    if is_element:
        _write_summary_data(coupling_summary_file, summary_data_f, summary_data_C)


def get_coupling_summary_file(save_path):
        coupling_summary_file = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbs_summary_cou.out"))

        coupling_summary_file.add_column_names(["NAME", "S",
                                                "F1001ABSPROP", "ERRF1001ABSPROP", "F1001REPROP", "ERRF1001REPROP", "F1001IMPROP", "ERRF1001IMPROP",
                                                "F1010ABSPROP", "ERRF1010ABSPROP", "F1010REPROP", "ERRF1010REPROP", "F1010IMPROP", "ERRF1010IMPROP",
                                                "PROPC11", "ERRPROPC11", "PROPC12", "ERRPROPC12", "PROPC21", "ERRPROPC21", "PROPC22", "ERRPROPC22",
                                                "MODEL_F1001RE", "MODEL_F1001IM", "MODEL_F1010RE", "MODEL_F1010IM",
                                                "MODEL_S"
                                                ])
        coupling_summary_file.add_column_datatypes(["%s", "%le",
                                                    "%le", "%le", "%le", "%le", "%le", "%le",
                                                    "%le", "%le", "%le", "%le", "%le", "%le",
                                                    "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le",
                                                    "%le", "%le", "%le", "%le",
                                                    "%le"])

        return coupling_summary_file


def _write_summary_data(coupling_summary_file, summary_data_f, summary_data_C):
    coupling_summary_file.add_table_row([summary_data_f[0], summary_data_f[1],
                                         summary_data_f[2], summary_data_f[3], summary_data_f[4], summary_data_f[5], summary_data_f[6], summary_data_f[7],
                                         summary_data_f[8], summary_data_f[9], summary_data_f[10], summary_data_f[11], summary_data_f[12], summary_data_f[13],
                                         summary_data_C[2], summary_data_C[3], summary_data_C[4], summary_data_C[5], summary_data_C[6], summary_data_C[7], summary_data_C[8], summary_data_C[9],
                                         summary_data_f[14], summary_data_f[15], summary_data_f[16], summary_data_f[17],
                                         summary_data_f[18]])


def _get_coupling_tfs_files(element_name, save_path, is_element):
    file_f_terms = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbscouple_" + element_name + ".out"))  # TODO: Probably we want a better name for these files
    file_C_terms = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbscouple2_" + element_name + ".out"))  # TODO: Probably we want a better name for these files

    if not is_element:
        file_f_terms.add_column_names(["NAME", "S",
                                       "F1001ABSMEAS", "ERRF1001ABSMEAS", "F1001REMEAS", "ERRF1001REMEAS", "F1001IMMEAS", "ERRF1001IMMEAS",
                                       "F1010ABSMEAS", "ERRF1010ABSMEAS", "F1010REMEAS", "ERRF1010REMEAS", "F1010IMMEAS", "ERRF1010IMMEAS",
                                       "F1001ABSCOR", "ERRF1001ABSCOR", "F1001RECOR", "ERRF1001RECOR", "F1001IMCOR", "ERRF1001IMCOR",
                                       "F1010ABSCOR", "ERRF1010ABSCOR", "F1010RECOR", "ERRF1010RECOR", "F1010IMCOR", "ERRF1010IMCOR",
                                       "F1001ABSBACK", "ERRF1001ABSBACK", "F1001REBACK", "ERRF1001REBACK", "F1001IMBACK", "ERRF1001IMBACK",
                                       "F1010ABSBACK", "ERRF1010ABSBACK", "F1010REBACK", "ERRF1010REBACK", "F1010IMBACK", "ERRF1010IMBACK",
                                       "F1001ABSBACKCOR", "ERRF1001ABSBACKCOR", "F1001REBACKCOR", "ERRF1001REBACKCOR", "F1001IMBACKCOR", "ERRF1001IMBACKCOR",
                                       "F1010ABSBACKCOR", "ERRF1010ABSBACKCOR", "F1010REBACKCOR", "ERRF1010REBACKCOR", "F1010IMBACKCOR", "ERRF1010IMBACKCOR",
                                       "MODEL_F1001RE", "MODEL_F1001IM", "MODEL_F1010RE", "MODEL_F1010IM", "MODEL_S"])
        file_f_terms.add_column_datatypes(["%s", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le"])

        file_C_terms.add_column_names(["NAME", "S",
                                       "PROPC11", "ERRPROPC11", "PROPC12", "ERRPROPC12", "PROPC21", "ERRPROPC21", "PROPC22", "ERRPROPC22",
                                       "CORC11", "ERRCORC11", "CORC12", "ERRCORC12", "CORC21", "ERRCORC21", "CORC22", "ERRCORC22",
                                       "BACKC11", "ERRBACKC11", "BACKC12", "ERRBACKC12", "BACKC21", "ERRBACKC21", "BACKC22", "ERRBACKC22",
                                       "BACKCORC11", "ERRBACKCORC11", "BACKCORC12", "ERRBACKCORC12", "BACKCORC21", "ERRBACKCORC21", "BACKCORC22", "ERRBACKCORC22",
                                       "MODEL_S"])
        file_C_terms.add_column_datatypes(["%s", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le"])
    else:
        file_f_terms.add_column_names(["NAME", "S",
                                       "F1001ABSPROP", "ERRF1001ABSPROP", "F1001REPROP", "ERRF1001REPROP", "F1001IMPROP", "ERRF1001IMPROP",
                                       "F1010ABSPROP", "ERRF1010ABSPROP", "F1010REPROP", "ERRF1010REPROP", "F1010IMPROP", "ERRF1010IMPROP",
                                       "MODEL_F1001RE", "MODEL_F1001IM", "MODEL_F1010RE", "MODEL_F1010IM", "MODEL_S"])
        file_f_terms.add_column_datatypes(["%s", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le", "%le", "%le", "%le", "%le"])

        file_C_terms.add_column_names(["NAME", "S",
                                       "PROPC11", "ERRPROPC11", "PROPC12", "ERRPROPC12", "PROPC21", "ERRPROPC21", "PROPC22", "ERRPROPC22",
                                       "MODEL_S"])
        file_C_terms.add_column_datatypes(["%s", "%le",
                                           "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le",
                                           "%le"])

    return file_f_terms, file_C_terms


def _write_f_terms(file_f_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor):
    summary_data = []
    first_bpm = bpms_list[0][1]
    last_bpm = bpms_list[-1][1]

    (f1001ab_ini, f1001_std_ini, p1001_ini, p1001_std_ini,
     f1001ab_end, f1001_std_end, p1001_end, p1001_std_end) = _get_start_end_f1001(measured_coupling, first_bpm, last_bpm)
    (f1010ab_ini, f1010_std_ini, p1010_ini, p1010_std_ini,
     f1010ab_end, f1010_std_end, p1010_end, p1010_std_end) = _get_start_end_f1010(measured_coupling, first_bpm, last_bpm)

    if input_model is not None:
        input_model.Cmatrix()

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        delta_phase_prop_x = model_propagation.MUX[model_propagation.indx[bpm_name]] % 1
        delta_phase_prop_y = model_propagation.MUY[model_propagation.indx[bpm_name]] % 1

        f1001_prop = model_propagation.f1001[model_propagation.indx[bpm_name]]
        err_f1001re_prop = _propagate_error_coupling_1001_re(f1001ab_ini, p1001_ini, delta_phase_prop_x, delta_phase_prop_y, f1001_std_ini, p1001_std_ini)
        err_f1001im_prop = _propagate_error_coupling_1001_im(f1001ab_ini, p1001_ini, delta_phase_prop_x, delta_phase_prop_y, f1001_std_ini, p1001_std_ini)
        err_f1001abs_prop = f1001_std_ini
        f1010_prop = model_propagation.f1010[model_propagation.indx[bpm_name]]
        err_f1010re_prop = _propagate_error_coupling_1010_re(f1010ab_ini, p1010_ini, delta_phase_prop_x, delta_phase_prop_y, f1010_std_ini, p1010_std_ini)
        err_f1010im_prop = _propagate_error_coupling_1010_im(f1010ab_ini, p1010_ini, delta_phase_prop_x, delta_phase_prop_y, f1010_std_ini, p1010_std_ini)
        err_f1010abs_prop = f1010_std_ini

        delta_phase_back_x = model_back_propagation.MUX[model_back_propagation.indx[bpm_name]] % 1
        delta_phase_back_y = model_back_propagation.MUY[model_back_propagation.indx[bpm_name]] % 1

        f1001_back = model_back_propagation.f1001[model_back_propagation.indx[bpm_name]]
        err_f1001re_back = _propagate_error_coupling_1001_re(f1001ab_end, p1001_end, delta_phase_back_x, delta_phase_back_y, f1001_std_end, p1001_std_end)
        err_f1001im_back = _propagate_error_coupling_1001_im(f1001ab_end, p1001_end, delta_phase_back_x, delta_phase_back_y, f1001_std_end, p1001_std_end)
        err_f1001abs_back = f1001_std_end
        f1010_back = model_back_propagation.f1010[model_back_propagation.indx[bpm_name]]
        err_f1010re_back = _propagate_error_coupling_1010_re(f1010ab_end, p1010_end, delta_phase_back_x, delta_phase_back_y, f1010_std_end, p1010_std_end)
        err_f1010im_back = _propagate_error_coupling_1010_im(f1010ab_end, p1010_end, delta_phase_back_x, delta_phase_back_y, f1010_std_end, p1010_std_end)
        err_f1010abs_back = f1010_std_end

        if not is_element:
            model_s = measured_coupling.S[measured_coupling.indx[bpm_name]]
            model_f1001r = measured_coupling.MDLF1001R[measured_coupling.indx[bpm_name]]
            model_f1001i = measured_coupling.MDLF1001I[measured_coupling.indx[bpm_name]]
            model_f1010r = measured_coupling.MDLF1010R[measured_coupling.indx[bpm_name]]
            model_f1010i = measured_coupling.MDLF1010I[measured_coupling.indx[bpm_name]]

            meas_abs_f1001 = measured_coupling.F1001W[measured_coupling.indx[bpm_name]]
            meas_q1001 = measured_coupling.Q1001[measured_coupling.indx[bpm_name]]
            meas_f1001r = measured_coupling.F1001R[measured_coupling.indx[bpm_name]]
            meas_f1001i = measured_coupling.F1001I[measured_coupling.indx[bpm_name]]
            err_meas_abs_f1001 = measured_coupling.FWSTD1[measured_coupling.indx[bpm_name]]
            err_meas_q1001 = measured_coupling.Q1001STD[measured_coupling.indx[bpm_name]]
            err_meas_f1001r = _propagate_error_coupling_1001_re(meas_abs_f1001, meas_q1001, 0., 0., err_meas_abs_f1001, err_meas_q1001)
            err_meas_f1001i = _propagate_error_coupling_1001_im(meas_abs_f1001, meas_q1001, 0., 0., err_meas_abs_f1001, err_meas_q1001)

            meas_abs_f1010 = measured_coupling.F1010W[measured_coupling.indx[bpm_name]]
            meas_q1010 = measured_coupling.Q1010[measured_coupling.indx[bpm_name]]
            meas_f1010r = measured_coupling.F1010R[measured_coupling.indx[bpm_name]]
            meas_f1010i = measured_coupling.F1010I[measured_coupling.indx[bpm_name]]
            err_meas_abs_f1010 = measured_coupling.FWSTD1[measured_coupling.indx[bpm_name]]
            err_meas_q1010 = measured_coupling.Q1010STD[measured_coupling.indx[bpm_name]]
            err_meas_f1010r = _propagate_error_coupling_1010_re(meas_abs_f1010, meas_q1010, 0., 0., err_meas_abs_f1010, err_meas_q1010)
            err_meas_f1010i = _propagate_error_coupling_1010_im(meas_abs_f1010, meas_q1010, 0., 0., err_meas_abs_f1010, err_meas_q1010)

            delta_phase_corr_x = model_cor.MUX[model_cor.indx[bpm_name]] % 1
            delta_phase_corr_y = model_cor.MUY[model_cor.indx[bpm_name]] % 1

            f1001_corr = model_cor.f1001[model_cor.indx[bpm_name]]
            err_f1001re_corr = _propagate_error_coupling_1001_re(f1001ab_ini, p1001_ini, delta_phase_corr_x, delta_phase_corr_y, f1001_std_ini, p1001_std_ini)
            err_f1001im_corr = _propagate_error_coupling_1001_im(f1001ab_ini, p1001_ini, delta_phase_corr_x, delta_phase_corr_y, f1001_std_ini, p1001_std_ini)
            err_f1001abs_corr = f1001_std_ini
            f1010_corr = model_cor.f1010[model_cor.indx[bpm_name]]
            err_f1010re_corr = _propagate_error_coupling_1010_re(f1010ab_ini, p1010_ini, delta_phase_corr_x, delta_phase_corr_y, f1010_std_ini, p1010_std_ini)
            err_f1010im_corr = _propagate_error_coupling_1010_im(f1010ab_ini, p1010_ini, delta_phase_corr_x, delta_phase_corr_y, f1010_std_ini, p1010_std_ini)
            err_f1010abs_corr = f1010_std_ini

            delta_phase_back_corr_x = model_back_cor.MUX[model_back_cor.indx[bpm_name]] % 1
            delta_phase_back_corr_y = model_back_cor.MUY[model_back_cor.indx[bpm_name]] % 1

            f1001_back_corr = model_back_cor.f1001[model_back_cor.indx[bpm_name]]
            err_f1001re_back_corr = _propagate_error_coupling_1001_re(f1001ab_end, p1001_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1001_std_end, p1001_std_end)
            err_f1001im_back_corr = _propagate_error_coupling_1001_im(f1001ab_end, p1001_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1001_std_end, p1001_std_end)
            err_f1001abs_back_corr = f1001_std_end
            f1010_back_corr = model_back_cor.f1010[model_back_cor.indx[bpm_name]]
            err_f1010re_back_corr = _propagate_error_coupling_1010_re(f1010ab_end, p1010_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1010_std_end, p1010_std_end)
            err_f1010im_back_corr = _propagate_error_coupling_1010_im(f1010ab_end, p1010_end, delta_phase_back_corr_x, delta_phase_back_corr_y, f1010_std_end, p1010_std_end)
            err_f1010abs_back_corr = f1010_std_end

            file_f_terms.add_table_row([bpm_name, bpm_s,
                                        meas_abs_f1001, err_meas_abs_f1001, meas_f1001r, err_meas_f1001r, meas_f1001i, err_meas_f1001i,
                                        meas_abs_f1010, err_meas_abs_f1010, meas_f1010r, err_meas_f1010r, meas_f1010i, err_meas_f1010i,
                                        abs(f1001_corr), err_f1001abs_corr, f1001_corr.real, err_f1001re_corr, f1001_corr.imag, err_f1001im_corr,
                                        abs(f1010_corr), err_f1010abs_corr, f1010_corr.real, err_f1010re_corr, f1010_corr.imag, err_f1010im_corr,
                                        abs(f1001_back), err_f1001abs_back, f1001_back.real, err_f1001re_back, f1001_back.imag, err_f1001im_back,
                                        abs(f1010_back), err_f1010abs_back, f1010_back.real, err_f1010re_back, f1010_back.imag, err_f1010im_back,
                                        abs(f1001_back_corr), err_f1001abs_back_corr, f1001_back_corr.real, err_f1001re_back_corr, f1001_back_corr.imag, err_f1001im_back_corr,
                                        abs(f1010_back_corr), err_f1010abs_back_corr, f1010_back_corr.real, err_f1010re_back_corr, f1010_back_corr.imag, err_f1010im_back_corr,
                                       model_f1001r, model_f1001i, model_f1010r, model_f1010i, model_s])
        else:
            model_s = input_model.S[input_model.indx[bpm_name]]

            model_f1001r = input_model.f1001[input_model.indx[bpm_name]].real
            model_f1001i = input_model.f1001[input_model.indx[bpm_name]].imag
            model_f1010r = input_model.f1010[input_model.indx[bpm_name]].real
            model_f1010i = input_model.f1010[input_model.indx[bpm_name]].imag

            average_f1001ab, final_f1001abs_err = weighted_average_for_SbS_elements(abs(f1001_prop), err_f1001abs_prop, abs(f1001_back), err_f1001abs_back)
            average_f1001re, final_f1001re_err = weighted_average_for_SbS_elements(f1001_prop.real, err_f1001re_prop, f1001_back.real, err_f1001re_back)
            average_f1001im, final_f1001im_err = weighted_average_for_SbS_elements(f1001_prop.imag, err_f1001im_prop, f1001_back.imag, err_f1001im_back)

            average_f1010ab, final_f1010abs_err = weighted_average_for_SbS_elements(abs(f1010_prop), err_f1010abs_prop, abs(f1010_back), err_f1010abs_back)
            average_f1010re, final_f1010re_err = weighted_average_for_SbS_elements(f1010_prop.real, err_f1010re_prop, f1010_back.real, err_f1010re_back)
            average_f1010im, final_f1010im_err = weighted_average_for_SbS_elements(f1010_prop.imag, err_f1010im_prop, f1010_back.imag, err_f1010im_back)

            file_f_terms.add_table_row([bpm_name, bpm_s,
                                        average_f1001ab, final_f1001abs_err, average_f1001re, final_f1001re_err, average_f1001im, final_f1001im_err,
                                        average_f1010ab, final_f1010abs_err, average_f1010re, final_f1010re_err, average_f1010im, final_f1010im_err,
                                        model_f1001r, model_f1001i, model_f1010r, model_f1010i, model_s])
            if bpm_name == element_name:
                summary_data = [bpm_name, bpm_s,
                                average_f1001ab, final_f1001abs_err, average_f1001re, final_f1001re_err, average_f1001im, final_f1001im_err,
                                average_f1010ab, final_f1010abs_err, average_f1010re, final_f1010re_err, average_f1010im, final_f1010im_err,
                                model_f1001r, model_f1001i, model_f1010r, model_f1010i, model_s]
    file_f_terms.write_to_file()
    return summary_data


def _write_C_terms(file_C_terms, element_name, is_element, bpms_list, measured_coupling, input_model, model_propagation, model_back_propagation, model_cor, model_back_cor):

    summary_data = []

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        prop_c11, prop_c12, prop_c21, prop_c22 = model_propagation.C[model_propagation.indx[bpm_name]]
        err_prop_c11, err_prop_c12, err_prop_c21, err_prop_c22 = [1e-8, 1e-8, 1e-8, 1e-8]  # TODO: propagate errors

        back_c11, back_c12, back_c21, back_c22 = model_back_propagation.C[model_back_propagation.indx[bpm_name]]
        err_back_c11, err_back_c12, err_back_c21, err_back_c22 = [1e-8, 1e-8, 1e-8, 1e-8]  # TODO: propagate errors

        if not is_element:
            model_s = measured_coupling.S[measured_coupling.indx[bpm_name]]
            cor_c11, cor_c12, cor_c21, cor_c22 = model_cor.C[model_cor.indx[bpm_name]]
            err_cor_c11, err_cor_c12, err_cor_c21, err_cor_c22 = [1e-8, 1e-8, 1e-8, 1e-8]  # TODO: propagate errors

            back_cor_c11, back_cor_c12, back_cor_c21, back_cor_c22 = model_back_cor.C[model_back_cor.indx[bpm_name]]
            err_back_cor_c11, err_back_cor_c12, err_back_cor_c21, err_back_cor_c22 = [1e-8, 1e-8, 1e-8, 1e-8]  # TODO: propagate errors

            file_C_terms.add_table_row([bpm_name, bpm_s,
                                        prop_c11, err_prop_c11, prop_c12, err_prop_c12, prop_c21, err_prop_c21, prop_c22, err_prop_c22,
                                        cor_c11, err_cor_c11, cor_c12, err_cor_c12, cor_c21, err_cor_c21, cor_c22, err_cor_c22,
                                        back_c11, err_back_c11, back_c12, err_back_c12, back_c21, err_back_c21, back_c22, err_back_c22,
                                        back_cor_c11, err_back_cor_c11, back_cor_c12, err_back_cor_c12, back_cor_c21, err_back_cor_c21, back_cor_c22, err_back_cor_c22,
                                        model_s])
        else:
            model_s = input_model.S[input_model.indx[bpm_name]]
            average_c11, final_err_c11 = weighted_average_for_SbS_elements(prop_c11, err_prop_c11, back_c11, err_back_c11)
            average_c12, final_err_c12 = weighted_average_for_SbS_elements(prop_c12, err_prop_c12, back_c12, err_back_c12)
            average_c21, final_err_c21 = weighted_average_for_SbS_elements(prop_c21, err_prop_c21, back_c21, err_back_c21)
            average_c21, final_err_c21 = weighted_average_for_SbS_elements(prop_c22, err_prop_c22, back_c22, err_back_c22)
            file_C_terms.add_table_row([bpm_name, bpm_s,
                                        average_c11, final_err_c11, average_c12, final_err_c12, average_c21, final_err_c21, average_c21, final_err_c21,
                                        model_s])
            if bpm_name == element_name:
                summary_data = [bpm_name, bpm_s,
                                average_c11, final_err_c11, average_c12, final_err_c12, average_c21, final_err_c21, average_c21, final_err_c21,
                                model_s]
    file_C_terms.write_to_file()
    return summary_data


def _get_start_end_f1001(measured_coupling, first_bpm, last_bpm):
    f1001ab_ini = measured_coupling.F1001W[measured_coupling.indx[first_bpm]]
    f1001_std_ini = measured_coupling.FWSTD1[measured_coupling.indx[first_bpm]]
    p1001_ini = 2 * np.pi * measured_coupling.Q1001[measured_coupling.indx[first_bpm]]
    p1001_std_ini = 2 * np.pi * measured_coupling.Q1001STD[measured_coupling.indx[first_bpm]]

    f1001ab_end = measured_coupling.F1001W[measured_coupling.indx[last_bpm]]
    f1001_std_end = measured_coupling.FWSTD1[measured_coupling.indx[last_bpm]]
    p1001_end = 2 * np.pi * measured_coupling.Q1001[measured_coupling.indx[last_bpm]]
    p1001_std_end = 2 * np.pi * measured_coupling.Q1001STD[measured_coupling.indx[last_bpm]]

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
    return math.sqrt((f1001_std_ini * np.cos(2 * np.pi * (p1001_ini - phasex + phasey)))**2 + (2 * np.pi * p1001_std_ini * f1001ab_ini * np.sin(2 * np.pi * (p1001_ini - phasex + phasey)))**2)  # @IgnorePep8


def _propagate_error_coupling_1001_im(f1001ab_ini, p1001_ini, phasex, phasey, f1001_std_ini, p1001_std_ini):
    return math.sqrt((f1001_std_ini * np.sin(2 * np.pi * (p1001_ini - phasex + phasey)))**2 + (2 * np.pi * p1001_std_ini * f1001ab_ini * np.cos(2 * np.pi * (p1001_ini - phasex + phasey)))**2)  # @IgnorePep8


def _propagate_error_coupling_1010_re(f1010ab_ini, p1010_ini, phasex, phasey, f1010_std_ini, p1010_std_ini):
    return math.sqrt((f1010_std_ini * np.cos(2 * np.pi * (p1010_ini - phasex - phasey)))**2 + (2 * np.pi * p1010_std_ini * f1010ab_ini * np.sin(2 * np.pi * (p1010_ini - phasex - phasey)))**2)  # @IgnorePep8


def _propagate_error_coupling_1010_im(f1010ab_ini, p1010_ini, phasex, phasey, f1010_std_ini, p1010_std_ini):
    return math.sqrt((f1010_std_ini * np.sin(2 * np.pi * (p1010_ini - phasex - phasey)))**2 + (2 * np.pi * p1010_std_ini * f1010ab_ini * np.cos(2 * np.pi * (p1010_ini - phasex - phasey)))**2)  # @IgnorePep8
