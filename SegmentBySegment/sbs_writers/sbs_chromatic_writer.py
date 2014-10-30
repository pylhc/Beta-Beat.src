from Utilities import tfs_file_writer
import os
import SegmentBySegment


def write_chromatic(element_name, is_element, measured_chromatic_wx, measured_chromatic_wy, input_model, propagated_models, save_path, chrom_summary_file):
    file_chromatic_wx, file_chromatic_wy = _get_chromatic_w_files(save_path, element_name, is_element)

    model_propagation = propagated_models.propagation
    model_back_propagation = propagated_models.back_propagation
    model_cor = propagated_models.corrected
    model_back_cor = propagated_models.corrected_back_propagation

    if not is_element:
        bpms_list = SegmentBySegment.intersect([model_cor, model_propagation, model_back_propagation, model_back_cor, input_model, measured_chromatic_wx])
    else:
        bpms_list = SegmentBySegment.intersect([model_cor, model_propagation, model_back_propagation, model_back_cor, input_model])

    summary_data_x = _write_chromatic_w_for_plane(file_chromatic_wx, "X",
                                                  element_name, bpms_list, measured_chromatic_wx,
                                                  input_model, model_propagation, model_cor, model_back_propagation, model_back_cor,
                                                  save_path, is_element)
    summary_data_y = _write_chromatic_w_for_plane(file_chromatic_wy, "Y",
                                                  element_name, bpms_list, measured_chromatic_wy,
                                                  input_model, model_propagation, model_cor, model_back_propagation, model_back_cor,
                                                  save_path, is_element)
    if is_element:
        _write_summary_data(chrom_summary_file, summary_data_x, summary_data_y)


def get_chrom_summary_file(save_path):
        chrom_summary_file = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbs_summary_chrom.out"))
        chrom_summary_file.add_column_names(["NAME", "S",
                                             "WPROPX", "ERRWPROPX", "PHIPROPX", "ERRPHIPROPX",
                                             "WPROPY", "ERRWPROPY", "PHIPROPY", "ERRPHIPROPY",
                                             "MODELWX", "MODELWY", "MODELPHIX", "MODELPHIY", "MODEL_S"])
        chrom_summary_file.add_column_datatypes(["%bpm_s", "%le",
                                                 "%le", "%le", "%le", "%le",
                                                 "%le", "%le", "%le", "%le",
                                                 "%le", "%le", "%le", "%le", "%le"])
        return chrom_summary_file


def _get_chromatic_w_files(save_path, element_name, is_element):
    file_chromatic_wx = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetax_" + element_name + ".out"))
    file_chromatic_wy = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetay_" + element_name + ".out"))

    if not is_element:
        file_chromatic_wx.add_column_names(["NAME", "S",
                                            "WPROPX", "ERRWPROPX", "PHIPROPX", "ERRPHIPROPX",
                                            "WCORX", "ERRWCORX", "PHICORX", "ERRPHICORX",
                                            "WBACKX", "ERRWBACKX", "PHIBACKX", "ERRPHIBACKX",
                                            "WBACKCORX", "ERRWBACKCORX", "PHIBACKCORX", "ERRPHIBACKCORX",
                                            "MODELWX", "MODELPHIX", "MODEL_S"])
        file_chromatic_wx.add_column_datatypes(["%bpm_s", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le"])

        file_chromatic_wy.add_column_names(["NAME", "S",
                                            "WPROPY", "ERRWPROPY", "PHIPROPY", "ERRPHIPROPY",
                                            "WCORY", "ERRWCORY", "PHICORY", "ERRPHICORY",
                                            "WBACKY", "ERRWBACKY", "PHIBACKY", "ERRPHIBACKY",
                                            "WBACKCORY", "ERRWBACKCORY", "PHIBACKCORY", "ERRPHIBACKCORY",
                                            "MODELWY", "MODELPHIY", "MODEL_S"])
        file_chromatic_wy.add_column_datatypes(["%bpm_s", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le", "%le",
                                                "%le", "%le", "%le"])
    else:
        file_chromatic_wx.add_column_names(["NAME", "S", "WPROPX", "ERRWPROPX", "PHIPROPX", "ERRPHIPROPX", "MODELWX", "MODELPHIX", "MODEL_S"])
        file_chromatic_wx.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

        file_chromatic_wy.add_column_names(["NAME", "S", "WPROPY", "ERRWPROPY", "PHIPROPY", "ERRPHIPROPY", "MODELWY", "MODELPHIY", "MODEL_S"])
        file_chromatic_wy.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    return file_chromatic_wx, file_chromatic_wy


def _write_chromatic_w_for_plane(file_chromatic, plane, element_name, bpms_list, measured_chromatic_w,
                                    input_model, model_propagation, model_cor, model_back_propagation, model_back_cor,
                                    save_path, is_element):
    summary_data = []
    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        model_s = input_model.S[input_model.indx[bpm_name]]
        model_w = getattr(input_model, "W" + plane)[input_model.indx[bpm_name]]
        model_phi = getattr(input_model, "PHI" + plane)[input_model.indx[bpm_name]]

        w_prop = getattr(model_propagation, "W" + plane)[model_propagation.indx[bpm_name]]
        phi_prop = getattr(model_propagation, "PHI" + plane)[model_propagation.indx[bpm_name]]
        err_w_prop = 1e-8  # TODO: Propagate
        err_phi_prop = 1e-8  # TODO: Propagate

        w_back = getattr(model_back_propagation, "W" + plane)[model_back_propagation.indx[bpm_name]]
        phi_back = getattr(model_back_propagation, "PHI" + plane)[model_back_propagation.indx[bpm_name]]
        err_w_back = 1e-8  # TODO: Propagate
        err_phi_back = 1e-8  # TODO: Propagate

        if not is_element:
            w_cor = getattr(model_cor, "W" + plane)[model_cor.indx[bpm_name]]
            phi_cor = getattr(model_cor, "PHI" + plane)[model_cor.indx[bpm_name]]
            err_w_cor = 1e-8  # TODO: Propagate
            err_phi_cor = 1e-8  # TODO: Propagate

            w_back_cor = getattr(model_back_cor, "W" + plane)[model_back_cor.indx[bpm_name]]
            phi_back_cor = getattr(model_back_cor, "PHI" + plane)[model_back_cor.indx[bpm_name]]
            err_w_back_cor = 1e-8  # TODO: Propagate
            err_phi_back_cor = 1e-8  # TODO: Propagate

            file_chromatic.add_table_row([bpm_name, bpm_s,
                                          w_prop, err_w_prop, phi_prop, err_phi_prop,
                                          w_cor, err_w_cor, phi_cor, err_phi_cor,
                                          w_back, err_w_back, phi_back, err_phi_back,
                                          w_back_cor, err_w_back_cor, phi_back_cor, err_phi_back_cor,
                                          model_w, model_phi, model_s])
        else:
            average_w, final_err_w = SegmentBySegment.weighted_average_for_SbS_elements(w_prop, err_w_prop, w_back, err_w_prop)
            average_phi, final_err_phi = SegmentBySegment.weighted_average_for_SbS_elements(phi_prop, err_phi_prop, phi_back, err_phi_prop)
            file_chromatic.add_table_row([bpm_name, bpm_s, average_w, final_err_w, average_phi, final_err_phi, model_w, model_phi, model_s])
            if element_name == bpm_name:
                summary_data = [bpm_name, bpm_s, average_w, final_err_w, average_phi, final_err_phi, model_w, model_phi, model_s]

    file_chromatic.write_to_file()
    return summary_data


def  _write_summary_data(chrom_summary_file, summary_data_x, summary_data_y):
    chrom_summary_file.add_table_row([summary_data_x[0], summary_data_x[1],
                                      summary_data_x[2], summary_data_x[3], summary_data_x[4], summary_data_x[5],
                                      summary_data_y[2], summary_data_y[3], summary_data_y[4], summary_data_y[5],
                                      summary_data_x[6], summary_data_y[6], summary_data_x[7], summary_data_y[7], summary_data_x[8]])
