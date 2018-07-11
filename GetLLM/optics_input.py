import argparse
import sys
from model import manager


def parse_args():
    optics_input = OpticsInput.init_from_options(*_get_optics_parser())
    python_path = "/afs/cern.ch/work/o/omc/anaconda/bin/python "
    return optics_input, (python_path + " ".join(sys.argv))


class OpticsInput(object):
    DEFAULTS = {
        "max_closed_orbit": 4.0,
        "coupling_method": 2,
        "orbit_unit": "mm",
        "range_of_bpms": 11,
        "bpm_combinations": 10,
        "beta_model_cut": 0.15,
        "beta_error_cut": 0.15,
        "nprocesses": 16,
        "no_averaged_tune": False,
        "end_lattice_phase": False,
        "union": False,
        "nonlinear": False,
        "three_bpm_method": False,
        "only_coupling": False
    }

    def __init__(self):
        self.files = None
        self.outputdir = None
        self.max_closed_orbit = OpticsInput.DEFAULTS["max_closed_orbit"]
        self.coupling_method = OpticsInput.DEFAULTS["coupling_method"]
        self.orbit_unit = OpticsInput.DEFAULTS["orbit_unit"]
        self.range_of_bpms = OpticsInput.DEFAULTS["range_of_bpms"]
        self.bpm_combinations = OpticsInput.DEFAULTS["bpm_combinations"]
        self.beta_model_cut = OpticsInput.DEFAULTS["beta_model_cut"]
        self.beta_error_cut = OpticsInput.DEFAULTS["beta_error_cut"]
        self.nprocesses = OpticsInput.DEFAULTS["nprocesses"]
        self.no_averaged_tune = OpticsInput.DEFAULTS["no_averaged_tune"]
        self.end_lattice_phase = OpticsInput.DEFAULTS["end_lattice_phase"]
        self.union = OpticsInput.DEFAULTS["union"]
        self.nonlinear = OpticsInput.DEFAULTS["nonlinear"]
        self.three_bpm_method = OpticsInput.DEFAULTS["three_bpm_method"]
        self.only_coupling = OpticsInput.DEFAULTS["only_coupling"]
        self.accelerator = None

    @staticmethod
    def init_from_options(options, accelerator):
        self = OpticsInput()
        self.files = options.files
        self.outputdir = options.outputdir
        self.max_closed_orbit = options.max_closed_orbit
        self.coupling_method = options.coupling_method
        self.orbit_unit = options.orbit_unit
        self.range_of_bpms = options.range_of_bpms
        self.bpm_combinations = options.bpm_combinations
        self.beta_model_cut = options.beta_model_cut
        self.beta_error_cut = options.beta_error_cut
        self.nprocesses = options.nprocesses
        self.no_averaged_tune = options.no_averaged_tune
        self.end_lattice_phase = options.end_lattice_phase
        self.union = options.union
        self.nonlinear = options.nonlinear
        self.three_bpm_method = options.three_bpm_method
        self.only_coupling = options.only_coupling
        self.accelerator = accelerator
        return self


def _get_optics_parser(start_args=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    parser.add_argument("--files", dest="files", help="Files from analysis, separated by comma")
    parser.add_argument("--outputdir", dest="outputdir", help="Output directory")
    parser.add_argument("--max_closed_orbit", dest="max_closed_orbit",
                        default=OpticsInput.DEFAULTS["max_closed_orbit"],
                        help="Maximal closed orbit for dispersion measurement in 'orbit_unit'")
    parser.add_argument("--coupling_method", dest="coupling_method", type=int, choices=(1, 2),
                        default=OpticsInput.DEFAULTS["coupling_method"],
                        help="Analysis option for coupling, 1 BPM or 2 BPMs")
    parser.add_argument("--orbit_unit", dest="orbit_unit", type=str,
                        choices=("um", "mm", "cm", "m"),
                        default=OpticsInput.DEFAULTS["orbit_unit"], help="Unit of orbit position.")
    parser.add_argument("--range_of_bpms", dest="range_of_bpms", choices=(3, 5, 7, 9, 11, 13, 15),
                        type=int, default=OpticsInput.DEFAULTS["range_of_bpms"],
                        help="Range of BPMs for beta from phase calculation")
    parser.add_argument("--bpm_combinations", dest="bpm_combinations", type=int,
                        default=OpticsInput.DEFAULTS["bpm_combinations"],
                        help="Number of different BPM combinations for beta-calculation")
    parser.add_argument("--beta_model_cut", dest="beta_model_cut", type=float,
                        default=OpticsInput.DEFAULTS["beta_model_cut"],
                        help="Set beta-beating threshold for action calculations")
    parser.add_argument("--beta_error_cut", dest="beta_error_cut", type=float,
                        default=OpticsInput.DEFAULTS["beta_error_cut"],
                        help="Set beta relative uncertainty threshold for action calculations")
    # parser.add_argument("--calibration", dest="calibration_dir_path", help="Path to the directory
    # where the calibration files (calibration_x.out, calibration_y.out) are stored.", default=CALIBRATION)
    parser.add_argument("--nprocesses", dest="nprocesses", type=int,
                        default=OpticsInput.DEFAULTS["nprocesses"],
                        help="""Sets the number of processes used. -1: take the number of CPUs 0: run serially >1: take the specified number. default = {0:d}""".format(
                            OpticsInput.DEFAULTS["nprocesses"]))
    parser.add_argument("--no_averaged_tune", dest="no_averaged_tune", action="store_true",
                        help="The phases and amplitudes are calculated for tune found per BPM.")
    parser.add_argument("--end_lattice_phase", dest="end_lattice_phase", action="store_true",
                        help="Compensate phase shifts at the end of the lattice by tunes")
    parser.add_argument("--union", dest="union", action="store_true",
                        help="The phase per BPM is calculated from at least 3 valid measurements.")
    parser.add_argument("--nonlinear", dest="nonlinear", action="store_true",
                        help="Run the RDT analysis")
    parser.add_argument("--three_bpm_method", dest="three_bpm_method", action="store_true",
                        help="Forces to use the 3 BPM method")  # TODO do we need this or we put --no_systematic_errors
    parser.add_argument("--only_coupling", dest="only_coupling", action="store_true",
                        help="Only coupling is calculated. ")

    options, acc_args = parser.parse_known_args(args=start_args)

    accelerator = manager.get_accel_instance(acc_args)

    return options, accelerator


"""Some leftovers
   getllm_d.onlycoupling = onlycoupling
    
    def __init__(self):
        self.errordefspath = ""
        self.parallel = False
        self.important_pairs = {}
        

    

"""
