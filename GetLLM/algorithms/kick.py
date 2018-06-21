from model.accelerators.accelerator import AccExcitationMode
from utils import tfs_pandas
from utils import logging_tools
import logging

def _calculate_kick(getllm_d, twiss_d, phase_d, beta_d, files_dict, bbthreshold, errthreshold):
    '''
    Fills the following TfsFiles:
     - getkick.out
     - getkickac.out

    :returns: dict string --> GetllmTfsFile -- The same instace of files_dict to indicate that the dict was extended
    '''
    accelerator = getllm_d.accelerator

    mad_twiss = accelerator.get_model_tfs()
    if accelerator.excitation != AccExcitationMode.FREE:
        mad_ac = accelerator.get_driven_tfs()

    LOGGER.info( "Calculating kick")
    files = [twiss_d.zero_dpp_x + twiss_d.non_zero_dpp_x, twiss_d.zero_dpp_y + twiss_d.non_zero_dpp_y]
    common_index = twiss_d.non_zero_dpp_commonbpms_x.index.intersection(
        twiss_d.non_zero_dpp_commonbpms_y.index.intersection(
            beta_d.x_phase.keys())).intersection(beta_d.y_phase.keys())

    meansqrt_2jx = {}
    meansqrt_2jy = {}
    bpmrejx = {}
    bpmrejy = {}

    try:
        [meansqrt_2jx, meansqrt_2jy, _, _, tunes, dpp, bpmrejx, bpmrejy] = getkick(
            files, mad_twiss, beta_d, bbthreshold, errthreshold)
    except IndexError:  # occurs if either no x or no y files exist
        return files_dict, [], []

    #mean_2j = mean{2J} and meansqrt_2j=mean{sqrt(2J)}

    tfs_file_model = files_dict['getkick.out']
    tfs_file_model.add_comment("Calculates the kick from the model beta function")
    column_names_list = ["DPP", "QX", "QXRMS", "QY", "QYRMS", "NATQX", "NATQXRMS", "NATQY", "NATQYRMS", "sqrt2JX", "sqrt2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JX", "2JXSTD", "2JY", "2JYSTD"]
    column_types_list = ["%le", "%le", "%le", "%le", "%le",     "%le",      "%le",    "%le",      "%le", "%le",      "%le",        "%le",       "%le",    "%le",   "%le",  "%le",    "%le"]
    tfs_file_model.add_column_names(column_names_list)
    tfs_file_model.add_column_datatypes(column_types_list)

    for i in range(0, len(dpp)):
        list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i], tunes[4][i], tunes[5][i],
                            tunes[6][i], tunes[7][i], meansqrt_2jx['model'][i][0], meansqrt_2jx['model'][i][1],
                            meansqrt_2jy['model'][i][0], meansqrt_2jy['model'][i][1], (meansqrt_2jx['model'][i][0]**2),
                            (2*meansqrt_2jx['model'][i][0]*meansqrt_2jx['model'][i][1]),
                            (meansqrt_2jy['model'][i][0]**2),
                            (2*meansqrt_2jy['model'][i][0]*meansqrt_2jy['model'][i][1])]
        tfs_file_model.add_table_row(list_row_entries)
        actions_x, actions_y = meansqrt_2jx['phase'], meansqrt_2jy['phase']

    tfs_file_phase = files_dict['getkickphase.out']
    tfs_file_phase.add_float_descriptor("Threshold_for_abs(beta_d-beta_m)/beta_m", bbthreshold)
    tfs_file_phase.add_float_descriptor("Threshold_for_uncert(beta_d)/beta_d", errthreshold)
    tfs_file_phase.add_float_descriptor("X_BPMs_Rejected", bpmrejx['phase'][len(dpp) - 1])
    tfs_file_phase.add_float_descriptor("Y_BPMs_Rejected", bpmrejy['phase'][len(dpp) - 1])
    tfs_file_phase.add_column_names(column_names_list)
    tfs_file_phase.add_column_datatypes(column_types_list)
    for i in range(0, len(dpp)):
        list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i], tunes[4][i], tunes[5][i],
                            tunes[6][i], tunes[7][i], meansqrt_2jx['phase'][i][0], meansqrt_2jx['phase'][i][1],
                            meansqrt_2jy['phase'][i][0], meansqrt_2jy['phase'][i][1], (meansqrt_2jx['model'][i][0]**2),
                            (2*meansqrt_2jx['model'][i][0]*meansqrt_2jx['model'][i][1]),
                            (meansqrt_2jy['model'][i][0]**2),
                            (2*meansqrt_2jy['model'][i][0]*meansqrt_2jy['model'][i][1])]
        tfs_file_phase.add_table_row(list_row_entries)

    if getllm_d.accelerator.excitation != AccExcitationMode.FREE:
        tfs_file = files_dict['getkickac.out']
        tfs_file.add_float_descriptor("RescalingFactor_for_X", beta_d.x_ratio_f)
        tfs_file.add_float_descriptor("RescalingFactor_for_Y", beta_d.y_ratio_f)
        tfs_file.add_column_names(column_names_list + ["sqrt2JXRES", "sqrt2JXSTDRES", "sqrt2JYRES", "sqrt2JYSTDRES", "2JXRES", "2JXSTDRES", "2JYRES", "2JYSTDRES"])
        tfs_file.add_column_datatypes(column_types_list + ["%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        [inv_jx, inv_jy, tunes, dpp] = algorithms.compensate_excitation.getkickac(
            mad_ac, files, phase_d.ac2bpmac_x, phase_d.ac2bpmac_y, getllm_d.accelerator.get_beam_direction(), getllm_d.lhc_phase)
        for i in range(0, len(dpp)):
            #TODO: in table will be the ratio without f(beta_d.x_ratio) used but rescaling factor is f version(beta_d.x_ratio_f). Check it (vimaier)
            list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i], tunes[4][i], tunes[5][i], tunes[6][i], tunes[7][i], inv_jx[i][0], inv_jx[i][1], inv_jy[i][0], inv_jy[i][1], (inv_jx[i][0] ** 2), (2 * inv_jx[i][0] * inv_jx[i][1]), (inv_jy[i][0] ** 2), (2 * inv_jy[i][0] * inv_jy[i][1]), (inv_jx[i][0] / math.sqrt(beta_d.x_ratio)), (inv_jx[i][1] / math.sqrt(beta_d.x_ratio)), (inv_jy[i][0] / math.sqrt(beta_d.y_ratio)), (inv_jy[i][1] / math.sqrt(beta_d.y_ratio)), (inv_jx[i][0] ** 2 / beta_d.x_ratio), (2 * inv_jx[i][0] * inv_jx[i][1] / beta_d.x_ratio), (inv_jy[i][0] ** 2 / beta_d.y_ratio), (2 * inv_jy[i][0] * inv_jy[i][1] / beta_d.y_ratio)]
            tfs_file.add_table_row(list_row_entries)
            actions_x, actions_y = inv_jx, inv_jx

    return files_dict, actions_x, actions_y
# END _calculate_kick -------------------------------------------------------------------------------


#---- finding kick
def getkick(files,mad_twiss,beta_d,bbthreshold,errthreshold):

    #One source for each beta-function source
    Sources = ['model','amp','phase']

    #To count number of BPM rej in each calculation
    bpmrejx={}
    bpmrejy={}

    #meansqrt_2j-> will be mean{sqrt{2Jx}} and mean_2J-> will be mean{2J} where mean is taken over all BPMs not rejected
    meansqrt_2jx={}
    meansqrt_2jy={}
    mean_2jx={}
    mean_2jy={}

    for source in Sources:
        bpmrejx[source] = []
        bpmrejy[source] = []
        meansqrt_2jx[source] = []
        meansqrt_2jy[source] = []
        mean_2jx[source] = []
        mean_2jy[source] = []

    tunex = []
    tuney = []
    tunexRMS = []
    tuneyRMS = []

    nat_tunex = []
    nat_tuney = []
    nat_tunexRMS = []
    nat_tuneyRMS = []

    dpp=[]
    #thresholds for rejecting, input by user
    bbthreshold = float(bbthreshold)
    errthreshold = float(errthreshold)

    for j in range(0,len(files[0])):
        tw_x = files[0][j]
        tw_y = files[1][j]
        #Loop uses gen_kick_calc to get action for each beta function source in each plane
        #Note that the result for source='amp' is not currently being used
        for source in Sources:
            meansqrt_2jx_temp, mean_2jx_temp, rejected_bpm_countx = gen_kick_calc([tw_x], mad_twiss, beta_d, source, 'H', bbthreshold, errthreshold)
            meansqrt_2jy_temp, mean_2jy_temp, rejected_bpm_county = gen_kick_calc([tw_y], mad_twiss, beta_d, source, 'V', bbthreshold, errthreshold)
            meansqrt_2jx[source].append(meansqrt_2jx_temp)
            meansqrt_2jy[source].append(meansqrt_2jy_temp)
            mean_2jx[source].append(mean_2jx_temp)
            mean_2jy[source].append(mean_2jy_temp)
            bpmrejx[source].append(rejected_bpm_countx)
            bpmrejy[source].append(rejected_bpm_county)


        dpp.append(getattr(tw_x, "DPP", 0.0))

        tunex.append(getattr(tw_x, "Q1", 0.0))
        tuney.append(getattr(tw_y, "Q2", 0.0))
        tunexRMS.append(getattr(tw_x, "Q1RMS", 0.0))
        tuneyRMS.append(getattr(tw_y, "Q2RMS", 0.0))

        nat_tunex.append(getattr(tw_x, "NATQ1", 0.0))
        nat_tuney.append(getattr(tw_y, "NATQ2", 0.0))
        nat_tunexRMS.append(getattr(tw_x, "NATQ1RMS", 0.0))
        nat_tuneyRMS.append(getattr(tw_y, "NATQ2RMS", 0.0))

    tune_values_list = [tunex, tunexRMS, tuney, tuneyRMS, nat_tunex, nat_tunexRMS, nat_tuney, nat_tuneyRMS]
    return [meansqrt_2jx, meansqrt_2jy, mean_2jx, mean_2jy, tune_values_list, dpp, bpmrejx, bpmrejy]


def gen_kick_calc(list_of_files,mad_twiss,beta_d,source, plane, bbthreshold,errthreshold):
    '''
    Gets action for a given beta function source and plane
     :Parameters:
        'list_of_files': list
            list of either linx or liny files input, produced by drive
        'mad_twiss': twiss
            The model file.
        'beta_d': BetaData
            Contains measured beta functions calculated in GetLLM
        'source': string
            From where the beta function was calcuated, either model, phase, or amp
        returns:  list, list, int
            The value and uncertainty for sqrt(2J) and 2J resp., along #BPM rej
    '''
    commonbpms = utils.bpm.intersect(list_of_files)
    if source == 'model':
        commonbpms = utils.bpm.model_intersect(commonbpms, mad_twiss)
    elif source == 'amp':
        if plane == 'H':
            commonbpms = utils.bpm.intersect_with_bpm_list(commonbpms, beta_d.x_amp.keys())
        elif plane == 'V':
            commonbpms = utils.bpm.intersect_with_bpm_list(commonbpms, beta_d.y_amp.keys())
    elif source == 'phase':
        if plane == 'H':
            commonbpms = utils.bpm.intersect_with_bpm_list(commonbpms, beta_d.x_phase.keys())
        elif plane == 'V':
            commonbpms = utils.bpm.intersect_with_bpm_list(commonbpms, beta_d.y_phase.keys())
    commonbpms = utils.bpm.get_list_of_tuples(commonbpms)

    meansqrt2j = []
    mean2j = []
    rejbpmcount = 0

    #Two different action calcs follow, one is by finding mean(sqrt(2j))-->meansqrt2j and the other is by finding mean(2j)-->mean2j
    for i in range(0, len(commonbpms)):
        bn1 = str.upper(commonbpms[i][1])
        #Gets beta function for use in action calculation- performs checks on its quality, rejects if above threshold
        if plane == 'H':
            if source == 'model':
                tembeta = mad_twiss.BETX[mad_twiss.indx[bn1]]
            elif source == 'amp':
                if beta_d.x_amp[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif abs((beta_d.x_amp[bn1][0] - mad_twiss.BETX[mad_twiss.indx[bn1]])/mad_twiss.BETX[mad_twiss.indx[bn1]]) > bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.x_amp[bn1][1]/beta_d.x_amp[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.x_amp[bn1][0]
            elif source == 'phase':
                if     beta_d.x_phase[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif abs((beta_d.x_phase[bn1][0] - mad_twiss.BETX[mad_twiss.indx[bn1]])/mad_twiss.BETX[mad_twiss.indx[bn1]]) > bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.x_phase[bn1][1]/beta_d.x_phase[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.x_phase[bn1][0]
        elif plane == 'V':
            if source == 'model':
                tembeta = mad_twiss.BETY[mad_twiss.indx[bn1]]
            elif source == 'amp':
                if     beta_d.y_amp[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif abs((beta_d.y_amp[bn1][0] - mad_twiss.BETY[mad_twiss.indx[bn1]])/mad_twiss.BETY[mad_twiss.indx[bn1]]) > bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.y_amp[bn1][1]/beta_d.y_amp[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.y_amp[bn1][0]
            elif source == 'phase':
                if     beta_d.y_phase[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif (abs(beta_d.y_phase[bn1][0] - mad_twiss.BETY[mad_twiss.indx[bn1]])/mad_twiss.BETY[mad_twiss.indx[bn1]])> bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.y_phase[bn1][1]/beta_d.y_phase[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.y_phase[bn1][0]

        meansqrt2j_i = 0.0
        mean2j_i = 0.0

        for tw_file in list_of_files:
            meansqrt2j_i += tw_file.PK2PK[tw_file.indx[bn1]] / 2.
            mean2j_i = (tw_file.PK2PK[tw_file.indx[bn1]] / 2)**2

        meansqrt2j_i = meansqrt2j_i / len(list_of_files)
        meansqrt2j.append(meansqrt2j_i / math.sqrt(tembeta))
        mean2j_i = mean2j_i / len(list_of_files)
        mean2j.append(mean2j_i / tembeta)

    if len(commonbpms) == rejbpmcount:
        print "Beta function calculated from", source, "in plane", plane, "is no good, no kick or action will be calculated from it."
        j_old, mean_2j = [0, 0], [0, 0]
        return j_old, mean_2j, rejbpmcount

    meansqrt2j = np.array(meansqrt2j)
    meansqrt2j_ave = np.average(meansqrt2j)
    meansqrt2j_std = math.sqrt(np.average(meansqrt2j*meansqrt2j)-meansqrt2j_ave**2+2.2e-16) #this is std dev, not RMS
    mean2j = np.array(mean2j)
    mean2j_ave = np.average(mean2j)
    mean2j_std = math.sqrt(np.average(mean2j*mean2j)-mean2j_ave**2+2.2e-16) #this is std dev, not RMS

    meansqrt_2j = [meansqrt2j_ave, meansqrt2j_std]
    mean_2j = [ mean2j_ave, mean2j_std]

    return meansqrt_2j, mean_2j, rejbpmcount
