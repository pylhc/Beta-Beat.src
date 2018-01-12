import __init__  # @UnusedImport
import os
from optparse import OptionParser
import sys
import multiprocessing
import time
from madx import madxrunner
from Utilities import iotools


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
NUM_PROCESSES = multiprocessing.cpu_count()  # Default number of processes to use in simulations


def _parse_args():
    parser = OptionParser()
    parser.add_option("-o", "--output",
                    help="Output directory for the results.",
                    metavar="output", default=CURRENT_PATH, dest="output_dir")
    parser.add_option("-p", "--processes",
                    help="Number of parallel processes to use in the simulation.",
                    metavar="numproc", default=NUM_PROCESSES, dest="num_processes")
    parser.add_option("-b", "--beam",
                    help="Beam to use, either LHCB1 or LHCB2.",
                    metavar="BEAM", dest="beam")
    parser.add_option("-e", "--energy",
                    help="The energy of the beam. It must be: 0.45TeV, 3.5TeV, 4TeV or 6.5TeV.",
                    metavar="ENERGY", dest="energy")
    options, _ = parser.parse_args()

    if options.beam is None:
        print >> sys.stderr, "Beam sequence must be defined, it must be LHCB1 or LHCB2."
        sys.exit(-1)
    beam = options.beam.upper().replace("LHC", "")
    if(beam not in ["B1", "B2"]):
        print >> sys.stderr, "Incorrect beam sequence, it must be LHCB1 or LHCB2."
        sys.exit(-1)
    if options.energy not in ["0.45TeV", "3.5TeV", "4TeV", "6.5TeV", "7TeV"]:
        print >> sys.stderr, "No valid energy specified, it must be one of: 0.45TeV, 3.5TeV, 4TeV, 6.5TeV or 7TeV"
        sys.exit(-1)
    if options.energy == "0.45TeV":
        options.energy = "inj"
    return options.output_dir, int(options.num_processes), beam, options.energy


def prepare_error_tables(output_dir, num_processes, beam, energy):
    error_tables_dir = os.path.abspath(os.path.join(output_dir, "error_tables_" + energy))
    iotools.create_dirs(error_tables_dir)

    try:
        os.symlink("/afs/cern.ch/eng/lhc/optics/V6.503", "db5")
    except(OSError):
        pass

    pool = multiprocessing.Pool(processes=num_processes)
    print "Preparing error files..."
    start_time = time.time()
    _parallel_prepare_error_files(error_tables_dir, beam, energy, pool)
    end_time = time.time()

    try:
        os.unlink("db5")
        os.unlink("ats")
    except(OSError):
        pass

    print "Done (" + str(end_time - start_time) + " seconds)\n"


def _parallel_prepare_error_files(error_tables_dir, beam, energy, pool):
    args = [(seed, error_tables_dir, beam, energy) for seed in range(1, 61)]
    tasks = pool.map_async(_prepare_single_error_file, args)
    tasks.wait()


def _prepare_single_error_file(seed_path_tuple):
    err_num = str((seed_path_tuple[0] % 60) + 1).zfill(4)
    error_tables_dir = seed_path_tuple[1]
    beam = seed_path_tuple[2]
    energy = seed_path_tuple[3]
    b2_errors_path = os.path.abspath(os.path.join(CURRENT_PATH, "..", "..", "MODEL", "LHCB", "b2_errors"))
    madx_job = ""
    with open(os.path.join(CURRENT_PATH, 'error_table.mask'), 'r') as infile:
        for line in infile:
            new_line = line
            new_line = new_line.replace("%ERR_NUM", err_num)
            new_line = new_line.replace("%RUN_DATA_PATH", error_tables_dir)
            new_line = new_line.replace("%BEAM", beam)
            new_line = new_line.replace("%ENERGY", energy)
            new_line = new_line.replace("%B2_ERRORS_PATH", b2_errors_path)
            new_line = new_line.replace("%NUM_BEAM", beam.replace("B", ""))
            madx_job += new_line + "\n"
    madxrunner.runForInputString(madx_job, stdout=open(os.devnull, "w"))


if __name__ == "__main__":
    _output_dir, _num_processes, _beam, _energy = _parse_args()
    prepare_error_tables(_output_dir, _num_processes, _beam, _energy)
