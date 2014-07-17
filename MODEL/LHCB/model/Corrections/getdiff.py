'''
Created on ???

@author: ???

@version: 1.0.1

getdiff.py will be executed after GetLLM.

Provide as first argument the path to the output files of GetLLM.

getdiff needs as input:
    - getbetax_free.out or getbetax.out
    - getbetay_free.out or getbetay.out
    - getcouple_free.out or getcouple.out
    - getDx.out
    - getDy.out
    - twiss_cor.dat
    - twiss_no.dat

getdiff produces the following output in the directory stated as argument:
    - bbx.out
    - bby.out
    - couple.out
    - dx.out
    - dy.out
    - phasex.out
    - phasey.out


Change history:
 - 1.0.1, vimaier, 04th June 2013:
     Added module description
     Reformatted into sections parse_args(), main() and main invocation
     Renamed variables
     Removed warnings from static analysis
'''


import os
import sys

import __init__  # @UnusedImport __init__ adds the path to Beta-Beat.src
import Python_Classes4MAD.metaclass
from numpy.ma.core import sqrt

#===================================================================================================
# parse_args()-function
#===================================================================================================


def parse_args():
    ''' Parses the sys.argv[1], checks for valid input and returns the path to src files  '''
    if len(sys.argv) == 2:
        path_to_src_files = sys.argv[1]
        if not os.path.isdir(path_to_src_files):
            print >> sys.stderr, "No valid directory:", path_to_src_files
            sys.exit(1)
        corrected_model_path = os.path.join(path_to_src_files, 'twiss_cor.dat')
        uncorrected_model_path = os.path.join(path_to_src_files, 'twiss_no.dat')
    elif len(sys.argv) == 4:
        path_to_src_files = sys.argv[1]
        if not os.path.isdir(path_to_src_files):
            print >> sys.stderr, "No valid directory:", path_to_src_files
            sys.exit(1)
        corrected_model_path = sys.argv[2]
        uncorrected_model_path = sys.argv[3]
    else:
        print >> sys.stderr, "Provide a path for a measurement file and optionally the path to the corrected and" \
                             "uncorrected models."
        sys.exit(1)

    return path_to_src_files, corrected_model_path, uncorrected_model_path

#===================================================================================================
# main()-function
#===================================================================================================


def main(path, corrected_model_path, uncorrected_model_path):
    '''
    :Parameters:
        'path': string
            Path to src files. Will also be used for output files.
    :Return: int
        0 if execution was successful otherwise !=0
    '''
    twiss_cor = Python_Classes4MAD.metaclass.twiss(corrected_model_path)
    twiss_no = Python_Classes4MAD.metaclass.twiss(uncorrected_model_path)
    twiss_cor.Cmatrix()
    twiss_no.Cmatrix()

    write_beta_diff_files(path, twiss_cor, twiss_no)

    write_coupling_diff_file(path, twiss_cor)

    write_dispersion_diff_files(path, twiss_cor, twiss_no)

    write_phase_diff_files(path, twiss_cor, twiss_no)

    write_chromatic_coupling_files(path, corrected_model_path)

    return 0


def write_beta_diff_files(path, twiss_cor, twiss_no):
    file_bbx = open(os.path.join(path, "bbx.out"), "w")
    file_bby = open(os.path.join(path, "bby.out"), "w")
    print >> file_bbx, "* NAME S MEA ERROR MODEL"
    print >> file_bbx, "$ %s %le %le %le %le"
    print >> file_bby, "* NAME S MEA ERROR MODEL"
    print >> file_bby, "$ %s %le %le %le %le"
    if os.path.exists(os.path.join(path, 'getbetax_free.out')):
        twiss_getbetax = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getbetax_free.out'))
    else:
        twiss_getbetax = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getbetax.out'))
    for i in range(len(twiss_getbetax.NAME)):
        bpm_name = twiss_getbetax.NAME[i]
        bpm_included = True
        try:
            check = twiss_cor.NAME[twiss_cor.indx[bpm_name]]  # @UnusedVariable
        except:
            print "No ", bpm_name
            bpm_included = False
        if bpm_included:
            j = twiss_cor.indx[bpm_name]
            t_x = twiss_getbetax  # Variable for abbreviation
            error_beta = sqrt(t_x.STDBETX[i] ** 2 + t_x.ERRBETX[i] ** 2) / t_x.BETXMDL[i]
            print >> file_bbx, bpm_name, t_x.S[i], (t_x.BETX[i] - t_x.BETXMDL[i]) / t_x.BETXMDL[i], error_beta, (twiss_cor.BETX[j] - twiss_no.BETX[j]) / twiss_no.BETX[j]

    if os.path.exists(os.path.join(path, 'getbetay_free.out')):
        twiss_getbetay = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getbetay_free.out'))
    else:
        twiss_getbetay = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getbetay.out'))
    for i in range(len(twiss_getbetay.NAME)):
        bpm_name = twiss_getbetay.NAME[i]
        bpm_included = True
        try:
            check = twiss_cor.NAME[twiss_cor.indx[bpm_name]]  # @UnusedVariable
        except:
            print "No ", bpm_name
            bpm_included = False
        if bpm_included:
            j = twiss_cor.indx[bpm_name]
            t_y = twiss_getbetay  # Variable for abbreviation
            error_beta = sqrt(t_y.STDBETY[i] ** 2 + t_y.ERRBETY[i] ** 2) / t_y.BETYMDL[i]
            print >> file_bby, bpm_name, t_y.S[i], (t_y.BETY[i] - t_y.BETYMDL[i]) / t_y.BETYMDL[i], error_beta, (twiss_cor.BETY[j] - twiss_no.BETY[j]) / twiss_no.BETY[j]

    file_bbx.close()
    file_bby.close()


def write_dispersion_diff_files(path, twiss_cor, twiss_no):
    try:
        twiss_getdx = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getDx.out'))
        file_dx = open(os.path.join(path, "dx.out"), "w")
        print >> file_dx, "* NAME S MEA ERROR MODEL"
        print >> file_dx, "$ %s %le %le %le %le"
        for i in range(len(twiss_getdx.NAME)):
            bpm_name = twiss_getdx.NAME[i]
            bpm_included = True
            try:
                check = twiss_cor.NAME[twiss_cor.indx[bpm_name]]  # @UnusedVariable
            except:
                print "No ", bpm_name
                bpm_included = False
            if bpm_included:
                j = twiss_cor.indx[bpm_name]
                print >> file_dx, bpm_name, twiss_getdx.S[i], twiss_getdx.DX[i] - twiss_getdx.DXMDL[i], twiss_getdx.STDDX[i], twiss_cor.DX[j] - twiss_no.DX[j]

        file_dx.close()
    except IOError:
        print "NO dispersion"
    except AttributeError:
        print "Empty table in getDx.out?! NO dispersion"
    try:
        twiss_getdy = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getDy.out'))
        file_dy = open(os.path.join(path, "dy.out"), "w")
        print >> file_dy, "* NAME S MEA ERROR MODEL"
        print >> file_dy, "$ %s %le %le %le %le"
        for i in range(len(twiss_getdy.NAME)):
            bpm_name = twiss_getdy.NAME[i]
            bpm_included = True
            try:
                check = twiss_cor.NAME[twiss_cor.indx[bpm_name]]  # @UnusedVariable
            except:
                print "No ", bpm_name
                bpm_included = False
            if bpm_included:
                j = twiss_cor.indx[bpm_name]
                print >> file_dy, bpm_name, twiss_getdy.S[i], twiss_getdy.DY[i] - twiss_getdy.DYMDL[i], twiss_getdy.STDDY[i], twiss_cor.DY[j] - twiss_no.DY[j]

        file_dy.close()
    except IOError:
        print "NO dispersion."
    except AttributeError:
        print "Empty table in getDy.out?! NO dispersion"


def write_coupling_diff_file(path, twiss_cor):
    file_couple = open(os.path.join(path, "couple.out"), "w")
    print >> file_couple, "* NAME S F1001re F1001im F1001e F1001re_m F1001im_m"
    print >> file_couple, "$ %s %le %le %le %le %le %le"
    if os.path.exists(os.path.join(path, 'getcouple_free.out')):
        twiss_getcouple = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getcouple_free.out'))
    else:
        twiss_getcouple = Python_Classes4MAD.metaclass.twiss(os.path.join(path, 'getcouple.out'))
    for i in range(len(twiss_getcouple.NAME)):
        bpm_name = twiss_getcouple.NAME[i]
        bpm_included = True
        try:
            check = twiss_cor.NAME[twiss_cor.indx[bpm_name]]  # @UnusedVariable
        except:
            print "No ", bpm_name
            bpm_included = False
        if bpm_included:
            j = twiss_cor.indx[bpm_name]
            print >> file_couple, bpm_name, twiss_getcouple.S[i], twiss_getcouple.F1001R[i], twiss_getcouple.F1001I[i], twiss_getcouple.FWSTD1[i], twiss_cor.f1001[j].real, twiss_cor.f1001[j].imag

    file_couple.close()


def write_phase_diff_files(path, twiss_cor, twiss_no):
    file_phase_x = open(os.path.join(path, "phasex.out"), "w")
    file_phase_y = open(os.path.join(path, "phasey.out"), "w")

    print >> file_phase_x, "* NAME S MEA ERROR MODEL DIFF DIFF_MDL"
    print >> file_phase_x, "$ %s %le %le %le %le %le %le"

    print >> file_phase_y, "* NAME S MEA ERROR MODEL DIFF DIFF_MDL"
    print >> file_phase_y, "$ %s %le %le %le %le %le %le"

    twiss_phase_x = _try_to_load_twiss(path, 'getphasex_free.out')
    twiss_phase_y = _try_to_load_twiss(path, 'getphasey_free.out')

    if not twiss_phase_x or not twiss_phase_y:
        twiss_phase_x = _try_to_load_twiss(path, 'getphasex.out')
        twiss_phase_y = _try_to_load_twiss(path, 'getphasey.out')

    if not twiss_phase_x or not twiss_phase_y:
        print "Cannot read phase files. NO phase."
        return

    common_bpm_names = _get_bpms_in_experiment_and_model(twiss_phase_x, twiss_cor)
    for i in range(len(common_bpm_names)):
        bpm_name = common_bpm_names[i]
        twiss_exp_bpm_indx = twiss_phase_x.indx[bpm_name]
        twiss_corr_bpm_indx = twiss_cor.indx[bpm_name]
        if i + 1 < len(common_bpm_names):
            next_bpm_name = common_bpm_names[i + 1]
            twiss_corr_next_bpm_indx = twiss_cor.indx[next_bpm_name]
            phase_corr = twiss_cor.MUX[twiss_corr_next_bpm_indx] - twiss_cor.MUX[twiss_corr_bpm_indx]
        else:  # last bpm, we have to take the first as next and subtract the total phase advance
            total_phs_adv = twiss_cor.Q1
            twiss_corr_first_bpm_indx = twiss_cor.indx[common_bpm_names[0]]
            phase_corr = twiss_cor.MUX[twiss_corr_first_bpm_indx] - twiss_cor.MUX[twiss_corr_bpm_indx] + total_phs_adv
        print >> file_phase_x, bpm_name,\
                    twiss_phase_x.S[twiss_exp_bpm_indx], twiss_phase_x.PHASEX[twiss_exp_bpm_indx],\
                    twiss_phase_x.STDPHX[twiss_exp_bpm_indx], phase_corr,\
                    twiss_phase_x.PHASEX[twiss_exp_bpm_indx] - twiss_phase_x.PHXMDL[twiss_exp_bpm_indx],\
                    phase_corr - twiss_phase_x.PHXMDL[twiss_exp_bpm_indx]

    common_bpm_names = _get_bpms_in_experiment_and_model(twiss_phase_y, twiss_cor)
    for i in range(len(common_bpm_names)):
        bpm_name = common_bpm_names[i]
        twiss_exp_bpm_indx = twiss_phase_y.indx[bpm_name]
        twiss_corr_bpm_indx = twiss_cor.indx[bpm_name]
        if i + 1 < len(common_bpm_names):
            next_bpm_name = common_bpm_names[i + 1]
            twiss_corr_next_bpm_indx = twiss_cor.indx[next_bpm_name]
            phase_corr = twiss_cor.MUY[twiss_corr_next_bpm_indx] - twiss_cor.MUY[twiss_corr_bpm_indx]
        else:  # last bpm, we have to take the first as next and subtract the total phase advance
            total_phs_adv = twiss_cor.Q2
            twiss_corr_first_bpm_indx = twiss_cor.indx[common_bpm_names[0]]
            phase_corr = twiss_cor.MUY[twiss_corr_first_bpm_indx] - twiss_cor.MUY[twiss_corr_bpm_indx] + total_phs_adv
        print >> file_phase_y, bpm_name, twiss_phase_y.S[twiss_exp_bpm_indx],\
                    twiss_phase_y.PHASEY[twiss_exp_bpm_indx], twiss_phase_y.STDPHY[twiss_exp_bpm_indx], phase_corr,\
                    twiss_phase_y.PHASEY[twiss_exp_bpm_indx] - twiss_phase_y.PHYMDL[twiss_exp_bpm_indx],\
                    phase_corr - twiss_phase_y.PHYMDL[twiss_exp_bpm_indx]

    file_phase_x.close()
    file_phase_y.close()


def write_chromatic_coupling_files(path, corrected_model_path):

    twiss_cor_plus = _try_to_load_twiss(os.path.split(corrected_model_path)[0], "twiss_cor_dpp.dat")
    twiss_cor_minus = _try_to_load_twiss(os.path.split(corrected_model_path)[0], "twiss_cor_dpm.dat")
    meas_chrom_coupling = _try_to_load_twiss(path, "chromcoupling.out")

    if not twiss_cor_plus or not twiss_cor_minus or not meas_chrom_coupling:
        print "Cannot read chromatic coupling files. NO chromatic coupling"
        return

    twiss_cor_plus.Cmatrix()
    twiss_cor_minus.Cmatrix()

    corr_cf1001r = []
    corr_cf1001i = []
    corr_cf1010r = []
    corr_cf1010i = []

    deltap = abs(twiss_cor_plus.DELTAP - twiss_cor_minus.DELTAP)

    for i in range(len(twiss_cor_plus.NAME)):
        ccoupling = (twiss_cor_plus.F1001R[i] - twiss_cor_minus.F1001R[i]) / deltap
        corr_cf1001r.append(ccoupling)

        ccoupling = (twiss_cor_plus.F1001I[i] - twiss_cor_minus.F1001I[i]) / deltap
        corr_cf1001i.append(ccoupling)

        ccoupling = (twiss_cor_plus.F1001R[i] - twiss_cor_minus.F1001R[i]) / deltap
        corr_cf1010r.append(ccoupling)

        ccoupling = (twiss_cor_plus.F1010I[i] - twiss_cor_minus.F1010I[i]) / deltap
        corr_cf1010i.append(ccoupling)

    file_chromatic_coupling = open(os.path.join(path, "chromatic_coupling.out"), "w")

    print >> file_chromatic_coupling, "* NAME S Cf1001r Cf1001rERR Cf1001i Cf1001iERR Cf1001r_model Cf1001i_model"
    print >> file_chromatic_coupling, "$ %s %le %le %le %le %le %le %le"

    common_bpm_names = _get_bpms_in_experiment_and_model(meas_chrom_coupling, twiss_cor_plus)
    for i in range(len(common_bpm_names)):
        bpm_name = common_bpm_names[i]
        twiss_exp_bpm_indx = meas_chrom_coupling.indx[bpm_name]
        twiss_corr_bpm_indx = twiss_cor_plus.indx[bpm_name]

        print >> file_chromatic_coupling, bpm_name, twiss_cor_plus.S[twiss_corr_bpm_indx],\
                    meas_chrom_coupling.Cf1001r[twiss_exp_bpm_indx], meas_chrom_coupling.Cf1001rERR[twiss_exp_bpm_indx],\
                    meas_chrom_coupling.Cf1001i[twiss_exp_bpm_indx], meas_chrom_coupling.Cf1001iERR[twiss_exp_bpm_indx],\
                    corr_cf1001r[twiss_corr_bpm_indx], corr_cf1001i[twiss_corr_bpm_indx]


def _try_to_load_twiss(twiss_path, twiss_file_name):
    twiss_full_path = os.path.join(twiss_path, twiss_file_name)
    if not os.path.isfile(twiss_full_path):
        print "Twiss file " + twiss_file_name + " doesn't exist"
        return False
    try:
        stderr = sys.stderr
        sys.stderr = sys.stdout  # Temporary redirect of standard error output to avoid annoying messages from twiss
        twiss_file = Python_Classes4MAD.metaclass.twiss(twiss_full_path)
        sys.stderr = stderr
    except:
        print "Can't read " + twiss_file_name + " twiss file"
        return False

    return twiss_file


def _get_bpms_in_experiment_and_model(experimental_twiss, model_twiss):
        common_bpms_list = []
        for bpm_name in model_twiss.NAME:
            if bpm_name.upper() in experimental_twiss.indx:
                common_bpms_list.append((experimental_twiss.S[experimental_twiss.indx[bpm_name]], bpm_name))
            else:
                print bpm_name, " in model but not in experiment"
        common_bpms_list.sort()
        return zip(*common_bpms_list)[1]

#===================================================================================================
# main invocation
#===================================================================================================


def _start():
    path_to_src_files, corrected_model_path, uncorrected_model_path = parse_args()

    print "Start getdiff.main..."
    return_value = main(path_to_src_files, corrected_model_path, uncorrected_model_path)
    print "getdiff.main finished with", return_value

    sys.exit(return_value)

if __name__ == "__main__":
    _start()
