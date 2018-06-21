import helper


def _calculate_orbit(getllm_d, twiss_d, tune_d, mad_twiss, files_dict):
    '''
    Calculates orbit and fills the following TfsFiles:
     - getCOx.out
     - getCOy.out
     - getCOx_dpp_' + str(k + 1) + '.out
     - getCOy_dpp_' + str(k + 1) + '.out

    :param _GetllmData getllm_d: accel is used(In-param, values will only be read)
    :param _TwissData twiss_d: Holds twiss instances of the src files. (In-param, values will only be read)
    :param _TuneData tune_d: Holds tunes and phase advances (In-param, values will only be read)

    :returns: (list, list, dict)
     - an list of dictionairies from horizontal computations
     - an list of dictionairies from vertical computations
     - the same dict as param files_dict to indicate that dict will be extended here.
    '''
    print 'Calculating orbit'
    list_of_co_x = []
    if twiss_d.has_zero_dpp_x():
        [cox, bpms] = calculate_orbit(mad_twiss, twiss_d.zero_dpp_x)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOx.out']
        tfs_file.add_string_descriptor("TABLE", 'ORBIT')
        tfs_file.add_string_descriptor("TYPE", 'ORBIT')
        # TODO: tfs_file.add_string_descriptor("SEQUENCE", getllm_d.accel)
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "X", "STDX", "XMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), cox[bn1][0], cox[bn1][1],
                                mad_twiss.loc[bn1, "X"], mad_twiss.loc[bn1, "MUX"]]
            tfs_file.add_table_row(list_row_entries)

        list_of_co_x.append(cox)
    list_of_co_y = []
    if twiss_d.has_zero_dpp_y():
        [coy, bpms] = calculate_orbit(mad_twiss, twiss_d.zero_dpp_y)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOy.out']
        tfs_file.add_string_descriptor("TABLE", 'ORBIT')
        tfs_file.add_string_descriptor("TYPE", 'ORBIT')
        #TODO: tfs_file.add_string_descriptor("SEQUENCE", getllm_d.accel)
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "Y", "STDY", "YMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), coy[bn1][0], coy[bn1][1],
                                mad_twiss.loc[bn1, "Y"], mad_twiss.loc[bn1, "MUY"]]
            tfs_file.add_table_row(list_row_entries)

        list_of_co_y.append(coy)
    #-------- Orbit for non-zero DPP
    if twiss_d.has_non_zero_dpp_x():
        k = 0
        for twiss_file in twiss_d.non_zero_dpp_x:
            list_with_single_twiss = []
            list_with_single_twiss.append(twiss_file)
            filename = 'getCOx_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = GetllmTfsFile(filename)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(twiss_file.filename)
            tfs_file.add_float_descriptor("DPP", float(twiss_file.DPP))
            tfs_file.add_float_descriptor("Q1", tune_d.q1)
            tfs_file.add_float_descriptor("Q2", tune_d.q2)
            [codpp, bpms] = calculate_orbit(mad_twiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "X", "STDX", "XMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = str.upper(bpms[i][1])
                bns1 = bpms[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), codpp[bn1][0], codpp[bn1][1],
                                    mad_twiss.loc[bn1, "X"], mad_twiss.loc[bn1, "MUX"]]
                tfs_file.add_table_row(list_row_entries)

            list_of_co_x.append(codpp)
            k += 1

    if twiss_d.has_non_zero_dpp_y():
        k = 0
        for twiss_file in twiss_d.non_zero_dpp_y:
            list_with_single_twiss = []
            list_with_single_twiss.append(twiss_file)
            filename = 'getCOy_dpp_' + str(k + 1) + '.out'
            files_dict[filename] = GetllmTfsFile(filename)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(twiss_file.filename)
            tfs_file.add_float_descriptor("DPP", float(twiss_file.DPP))
            tfs_file.add_float_descriptor("Q1", tune_d.q1)
            tfs_file.add_float_descriptor("Q2", tune_d.q2)
            [codpp, bpms] = calculate_orbit(mad_twiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "Y", "STDY", "YMDL", "MUYMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpms)):
                bn1 = str.upper(bpms[i][1])
                bns1 = bpms[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), codpp[bn1][0], codpp[bn1][1],
                                    mad_twiss.loc[bn1, "Y"], mad_twiss.loc[bn1, "MUY"]]
                tfs_file.add_table_row(list_row_entries)

            list_of_co_y.append(codpp)
            k += 1

    return list_of_co_x, list_of_co_y, files_dict
# END _calculate_orbit ------------------------------------------------------------------------------


def calculate_orbit(mad_twiss, list_of_files):

    commonbpms=utils.bpm.intersect(list_of_files)
    commonbpms=utils.bpm.model_intersect(commonbpms, mad_twiss)
    commonbpms = utils.bpm.get_list_of_tuples(commonbpms)
    co={} # Disctionary for output
    for i in range(0,len(commonbpms)):
        bpm_name=str.upper(commonbpms[i][1])
        coi=0.0
        coi2=0.0

        for tw_file in list_of_files:
            coi=coi + tw_file.CO[tw_file.indx[bpm_name]]
            coi2=coi2 + tw_file.CO[tw_file.indx[bpm_name]]**2

        coi = coi/len(list_of_files)
        corms = math.sqrt(coi2/len(list_of_files)-coi**2+2.2e-16)
        co[bpm_name] = [coi,corms]

    return [co, commonbpms]
