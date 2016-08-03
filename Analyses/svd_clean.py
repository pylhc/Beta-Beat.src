"""
.. module: Analyses.svd_clean

Created on 30 Dec 2003

Options::

use -h

Usage cmd line::

    python svd_clean.py -f file [other options]
    python svd_clean.py --file=Beam1@Turn@2012_03_15@04_19_55_001_0.sdds
    python svd_clean.py --file=Beam1@Turn@2012_03_15@04_19_55_001_0.sdds --turn=1 --p=0.00001 --sumsquare=0.925 --sing_val=57
    python svd_clean.py --file=Beam1@Turn@2012_03_15@04_19_55_001_0.sdds --output=cleaned.sdds

Usage in another Python module::

    import Analyses.svd_clean
    ...
    Analyses.svd_clean.clean_sdds_file("my_dirty_sdds_ascii_file.sdds")

.. moduleauthor:: Ram Calaga, Thomas Bach
"""

import sys
import os
import time
import argparse
import numpy
from numpy import dot as matrixmultiply
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD/")
import __init__  # @UnusedImport init will include paths
from metaclass import twiss
from datetime import datetime
from __builtin__ import max

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
BAD_BPMS_DIR = os.path.join(CURRENT_PATH, "bad_bpms_files")

# internal options
PRINT_TIMES = False  # If True, (more) execution times will be printed
PRINT_DEBUG = False  # If True, internal debug information will be printed

# internal constants
PLANE_X = "0"
PLANE_Y = "1"

# piotr: for data 20160408 01:45+ (40cm) B1
LIST_OF_KNOWN_BAD_BPMS = ["BPM.17L8.B1", "BPM.16L8.B1", "BPM.8R8.B1", "BPM.9R8.B1", # H B1 (big ones)
                          "BPM.26L8.B1","BPM.24L8.B1", "BPM.10R6.B1","BPM.8R6.B1","BPM.32R1.B1","BPM.33R1.B1", # V B1 (big ones)
                           "BPM.12R2.B1","BPM.13R2.B1","BPM.15R6.B1","BPM.16R6.B1","BPM.19L7.B1","BPM.18L7.B1", # H B1 (small ones)
                           "BPM.21R7.B1","BPM.22R7.B1","BPM.20R8.B1","BPM.21R8.B1","BPM.19L2.B1","BPM.18L2.B1",  # H B1 (small ones)
                           "BPMR.7L5.B1","BPM.6L5.B1","BPM.8L1.B1","BPM.6L1.B1"]

LIST_OF_WRONG_POLARITY_BPMS_BOTH_PLANES = []

RESYNC_FROM = datetime(2016, 4, 1)  # 1st of April of 2016 -> Date around when the BPM synchronization was changed.
LIST_OF_OUT_OF_SYNC_BPMS_B1 = ["BPM.33L2.B1", "BPM.32L2.B1", "BPM.31L2.B1", "BPM.30L2.B1",
                               "BPM.29L2.B1", "BPM.28L2.B1", "BPM.27L2.B1", "BPM.26L2.B1", "BPM.25L2.B1", "BPM.24L2.B1", "BPM.23L2.B1", "BPM.22L2.B1", "BPM.21L2.B1", "BPM.20L2.B1",
                               "BPM.19L2.B1", "BPM.18L2.B1", "BPM.17L2.B1", "BPM.16L2.B1", "BPM.15L2.B1", "BPM.14L2.B1", "BPM.13L2.B1", "BPM.12L2.B1", "BPM.11L2.B1", "BPM.10L2.B1",
                               "BPM.9L2.B1", "BPM.8L2.B1", "BPM.7L2.B1", "BPMR.6L2.B1", "BPMYB.5L2.B1", "BPMYB.4L2.B1", "BPMWI.4L2.B1", "BPMS.2L2.B1", "BPMSW.1L2.B1"]
LIST_OF_OUT_OF_SYNC_BPMS_B2 = ["BPM.34R8.B2", "BPM.33R8.B2", "BPM.32R8.B2", "BPM.31R8.B2", "BPM.30R8.B2",
                               "BPM.29R8.B2", "BPM.28R8.B2", "BPM.27R8.B2", "BPM.26R8.B2", "BPM.25R8.B2", "BPM.24R8.B2", "BPM.23R8.B2", "BPM.22R8.B2", "BPM.21R8.B2", "BPM.20R8.B2",
                               "BPM.19R8.B2", "BPM.18R8.B2", "BPM.17R8.B2", "BPM.16R8.B2", "BPM.15R8.B2", "BPM.14R8.B2", "BPM.13R8.B2", "BPM.12R8.B2", "BPM.11R8.B2", "BPM.10R8.B2",
                               "BPM.9R8.B2", "BPM.8R8.B2", "BPM_A.7R8.B2", "BPMR.6R8.B2", "BPMYB.5R8.B2", "BPMYB.4R8.B2", "BPMWI.4R8.B2", "BPMS.2R8.B2", "BPMSW.1R8.B2"]

#dev hints:
# print ">> Time for removeBadBpms: {}s".format(time.time() - time_start)  # does not work with python <= 2.6
# print ">> Time for removeBadBpms: {0}s".format(time.time() - time_start) # works with python 2.6
# -tbach


#===================================================================================================
# _parse_args()-function
#===================================================================================================

def _parse_args():
    """ Parses the arguments and returns args """
    parser = argparse.ArgumentParser(description="""Takes a SDDS ASCII file and cleans it.\n 1) remove BPM with bad data
    (constant, spikes, zero-values)\n 2) reduce noise with a SVD.""")

    parser.add_argument("-f", "--file", "--inputfile", 
        help="File to clean", 
        required=True, dest="inputfile")
    parser.add_argument("-n", "--newfile", "-o", "--outputfile", 
        help="File name for output file. Default (or empty) is override current file", 
        dest="outputfile")
    parser.add_argument("-t", "--turn", "--startturn",
        help="Turn number to start (first is 1, not 0). Default is first turn: %(default)s",
        default="1", dest="startturn", type=int)
    parser.add_argument("-m", "--maxturns",
        help="""Maximum number of turns to be analysed. Default is a number that is lower than the maximum 
             which can be handled by drive: %(default)s""",
        default="9500", dest="maxturns", type=int)
    parser.add_argument("-v", "--sing_val", "--singular_values_amount_to_keep",
        help="""Keep this amount of singular values in decreasing order, rest will be cut (set to 0). 
             Default is a large number: %(default)s""",
        default="100000", dest="singular_values_amount_to_keep", type=int)
    parser.add_argument("-p", "--p", "--pk-2-pk", "--min_peak_to_peak", 
        help="""Peak to peak amplitude cut. This removes BPMs where 
             abs(max(turn values) - min(turn values)) <= threshold. Default: %(default)s""",
        default="0.00000001", dest="min_peak_to_peak", type=float)
    parser.add_argument("-c", "--max-peak-cut", "--max-peak",
        help="""Maximum peak tolerance in mm. This removes BPMs where the maximum measured oscillation > threshold.
             Default: %(default)s""",
        default="20", dest="max_peak", type=float)
    parser.add_argument("-s", "--sumsquare", "--single_svd_bpm_threshold",
        help="""Threshold for single BPM dominating a mode. Should be > 0.9 for LHC. Default: %(default)s""",
        default="0.925", dest="single_svd_bpm_threshold", type=float)
    parser.add_argument("-r", "--resync",
        help="""Fix the synchronization of the BPMs if the file timestamp is present and is after May of 2016.""",
        dest="resync", action="store_true")
    parser.add_argument("-d", "--subtract",
        help="""Subtracts the closed orbit.""",
        dest="subtract_co", action="store_true")
    parser.add_argument("--use_test_mode",
        help="""Set testing mode, this prevents date writing. Useful for automated tests, because otherwise each run
        will have a new date, therefore a diff will fail""",
        action="store_true", dest="use_test_mode")
    parser.add_argument("--print_times",
        help="""Prints (more verbose) execution times""",
        action="store_true", dest="print_times")
    
    args = parser.parse_args()
    print str(args).replace("Namespace(", "Arguments:\n").replace(", ", "\n") + "\n"
    
#     sys.exit()
    if args.print_times:
        global PRINT_TIMES
        PRINT_TIMES = True
        
    return args


#===================================================================================================
# main-function
#===================================================================================================
def clean_sdds_file(options):
    """
    Invokes cleaning of given SDDS file

    For more details, see @_parse_args implementation or call this file with -h
    """

    _InputData.static_init(
       inputfile=options.inputfile,
       outputfile=options.outputfile,
       startturn_human=options.startturn,
       maxturns_human=options.maxturns,
       singular_values_amount_to_keep=options.singular_values_amount_to_keep,
       min_peak_to_peak=options.min_peak_to_peak,
       max_peak_cut=options.max_peak,
       single_svd_bpm_threshold=options.single_svd_bpm_threshold,
       resync=options.resync,
       subtract_co=options.subtract_co,
       use_test_mode=options.use_test_mode
    )

    start_time = time.time()
    _SvdHandler()
    print "Global Time: {0}s".format(time.time() - start_time)
    return 0


class _InputData(object):
    """This class holds all input variables for svd clean """

    @staticmethod
    def static_init(inputfile, outputfile, startturn_human, maxturns_human, singular_values_amount_to_keep, 
                    min_peak_to_peak, max_peak_cut, single_svd_bpm_threshold, resync, subtract_co, use_test_mode):
        _InputData.inputfile = inputfile
        print "inputfile:  " + str(inputfile)
        _InputData.outputfile = outputfile
        if not _InputData.outputfile or _InputData.outputfile == "":
            #_InputData.outputfile = str(inputfile) + "_svd_cleaned"
            _InputData.outputfile = _InputData.inputfile
            print "Warning: no outputfile given, output=input (*overwrite*)"
        print "outputfile: " + str(_InputData.outputfile)
        _InputData.startturn_human = startturn_human
        _InputData.startturn = startturn_human - 1
        _InputData.maxturns_human = maxturns_human
        _InputData.maxturns = maxturns_human - 1
        _InputData.singular_values_amount_to_keep = singular_values_amount_to_keep
        _InputData.min_peak_to_peak = min_peak_to_peak
        _InputData.max_peak_cut = max_peak_cut
        _InputData.single_svd_bpm_threshold = single_svd_bpm_threshold
        _InputData.resync = resync
        _InputData.subtract_co = subtract_co
        _InputData.use_test_mode = use_test_mode

    def __init__(self):
        raise NotImplementedError("static class _InputData cannot be instantiated")


class _SddsFile(object):
    """This class represents a SDDS file. It reads an existing SDDS file and creates & writes the new one"""
    
    def __init__(self, path_to_sddsfile):
        self.path_to_sddsfile = path_to_sddsfile
        self.bad_bpm_file = _BadBpmFile(path_to_sddsfile)
        self.parsed = False
        self.header = ""
        self.number_of_turns = 0
        self.dictionary_plane_to_bpms = {PLANE_X: _Bpms(), PLANE_Y: _Bpms()}

    def init(self):
        """Parses the current SDDS file and sets all member variables"""
        if (self.parsed):
            return

        try:
            swapped_bpms = twiss(os.path.join(BAD_BPMS_DIR, "SwappedBPMs_2016-06-20.tfs"))
        except IOError:
            print >> sys.stderr, "Can't read swapped BPMs file."
            swapped_bpms = None
        except ValueError:
            swapped_bpms = None

        last_number_of_turns = 0
        detected_number_of_turns = 0
        
        bpm_with_exact_zero_counter = 0
        known_bad_bpms_counter = 0
        flatbpm_counter = 0
        bpm_with_spike_counter = 0
        do_resync = False
        
        time_start = time.time()
        print "Extracting data from file " + str(self.path_to_sddsfile)

        
        with open(self.path_to_sddsfile, "r") as filesdds:
            for line in filesdds:  # Iterator over all lines -tbach
                # example: 0 BPMSY.4L1.B1 23461.055 <fp values whitespace separated>  -tbach
                
                if line.startswith("#"):  # we have a comment line -tbach
                    if line.startswith("#acquisition stamp :") and _InputData.resync:
                        raw_acq_stamp = line.replace("#acquisition stamp :", "")
                        try:
                            acq_stamp = int(float(raw_acq_stamp) / 1.0e9)  # The timestamp comes in nanoseconds...
                            acq_date = datetime.fromtimestamp(acq_stamp)
                            print "Acquisition date: ", str(acq_date)
                        except ValueError:
                            print >> sys.stderr, "Cannot parse timestamp, will not resynchronize the BPMs."
                        if acq_date is not None:
                            do_resync = RESYNC_FROM < acq_date
                    self.header += line
                    continue
                
                # detects if a turn has a value of exact zero and removes it since this was a workaround
                # for the large spikes and is still unwanted. Could possible remove a good BPM but unlikely -alangner
                if " 0.0 " in line: # example: 0 BPMSY.4L1.B1 23461.055 0.0 0.0 -tbach
                    reason_for_badbpm = "Found an exact zero"
                    badbpm = _BadBpm(line, reason_for_badbpm)
                    self.bad_bpm_file.add_badbpm(badbpm)
                    bpm_with_exact_zero_counter += 1
                    continue
                
                # from here, we have a data line -tbach
                list_splitted_line_values = line.split()
                number_of_columns = len(list_splitted_line_values)
                if (number_of_columns < 3):  # this is not a valid data line -tbach
                    print "not a valid data line: " + str(line)
                    continue
                
                plane = list_splitted_line_values[0]
                bpm_name = list_splitted_line_values[1]
                location = float(list_splitted_line_values[2])

                if swapped_bpms is not None:
                    if bpm_name in swapped_bpms.NAME:
                        if plane == swapped_bpms.PLANE[swapped_bpms.indx[bpm_name]]:
                            plane = swapped_bpms.PLANE2[swapped_bpms.indx[bpm_name]]
                            location = float(swapped_bpms.S2[swapped_bpms.indx[bpm_name]])
                            bpm_name = swapped_bpms.NAME2[swapped_bpms.indx[bpm_name]]

                if bpm_name in LIST_OF_KNOWN_BAD_BPMS: # Remove known bad BPMs
                    reason_for_badbpm = "removed (known) bad bpm: " + str(bpm_name)
                    #print reason_for_badbpm
                    badbpm = _BadBpm(line, reason_for_badbpm)
                    self.bad_bpm_file.add_badbpm(badbpm)
                    known_bad_bpms_counter += 1
                    continue 
    
                # first 3 entries metadata, rest should be turn data -tbach
                detected_number_of_turns = number_of_columns - 3
    
                if last_number_of_turns > 0:
                    if last_number_of_turns != detected_number_of_turns:
                        print """(plane {0}) BPM has a different number of turns then previous BPM. 
                                Current BPM: {1}, previous turns: {2}, current turns: {3}""" \
                                .format(plane, bpm_name, last_number_of_turns, detected_number_of_turns)
                        sys.exit(1)
                # the else part will happen only once, if last_number_of_turns <= 0
                # (or multiple times if we do not have any turn data) -tbach
                else:
                    if (_InputData.startturn_human > detected_number_of_turns):
                        print "startturn > detected_number_of_turns. startturn: {0}, detected_number_of_turns: {1}" \
                                .format(_InputData.startturn_human, detected_number_of_turns)
                        sys.exit(1)
                last_number_of_turns = detected_number_of_turns
    
                bpm_name_location_plane = (bpm_name, location, plane)
                # be careful with the bpm_name_location_plane format,
                # it is important for other parts of the program
                # -- tbach
    
                self.number_of_turns = min(_InputData.maxturns, detected_number_of_turns) - _InputData.startturn
                ndarray_line_data = numpy.array(\
                    list_splitted_line_values[3 + _InputData.startturn:3 + _InputData.startturn + self.number_of_turns], \
                    dtype=numpy.float64)
                
                if _InputData.subtract_co:                
                    ndarray_line_data = ndarray_line_data - numpy.average(ndarray_line_data)    #Removing the average orbit
                
                if bpm_name in LIST_OF_WRONG_POLARITY_BPMS_BOTH_PLANES:
                    ndarray_line_data = -1. * ndarray_line_data

                if do_resync:
                    if bpm_name in LIST_OF_OUT_OF_SYNC_BPMS_B1 or bpm_name in LIST_OF_OUT_OF_SYNC_BPMS_B2:
                        ndarray_line_data = numpy.roll(ndarray_line_data, -1)
                    if len(LIST_OF_OUT_OF_SYNC_BPMS_B1) != 0 or len(LIST_OF_OUT_OF_SYNC_BPMS_B2) != 0:
                        ndarray_line_data = numpy.delete(ndarray_line_data, len(ndarray_line_data) - 1)

                ndarray_line_data_max = numpy.max(ndarray_line_data)
                ndarray_line_data_min = numpy.min(ndarray_line_data)
                
                 # this block handles BPMs with the same values for all turns -tbach
                peak_to_peak_difference = numpy.abs(ndarray_line_data_max - ndarray_line_data_min)
                if peak_to_peak_difference <= _InputData.min_peak_to_peak:  # then do not use this BPM -tbach
                    reason_for_badbpm = "Flat BPM, the difference between min/max is smaller than " + str(_InputData.min_peak_to_peak)
                    badbpm = _BadBpm(bpm_name_location_plane, ndarray_line_data, reason_for_badbpm)
                    self.bad_bpm_file.add_badbpm(badbpm)
                    flatbpm_counter += 1
                    continue
                
                # detects the turn numbers of all occurrences of a spike > _InputData.max_peak_cut
                if max(abs(ndarray_line_data_min), abs(ndarray_line_data_max)) > _InputData.max_peak_cut:
                    reason_for_badbpm = "Found a spike: abs({0}) > {1}mm" \
                        .format(ndarray_line_data_max, max(abs(ndarray_line_data_min), abs(ndarray_line_data_max)))
                    badbpm = _BadBpm(bpm_name_location_plane, ndarray_line_data, reason_for_badbpm)
                    self.bad_bpm_file.add_badbpm(badbpm)
                    bpm_with_spike_counter += 1
                    continue
                
                self.dictionary_plane_to_bpms[plane].bpm_data.append(ndarray_line_data)
                self.dictionary_plane_to_bpms[plane].bpm_name_location_plane.append(bpm_name_location_plane)

        self.parsed = True

        print ""
        if bpm_with_exact_zero_counter > 0:
            print "Exact zeros detected. BPMs removed: " + str(bpm_with_exact_zero_counter)
        if known_bad_bpms_counter > 0:
            print "Known bad BPMs removed: " + str(known_bad_bpms_counter)
        if flatbpm_counter > 0:
            print "Flat BPMS detected (diff min/max <= {0}. BPMs removed: {1}".format(_InputData.min_peak_to_peak, flatbpm_counter)
        if bpm_with_spike_counter > 0:
            print "Spikes > {0}mm detected. BPMs removed: {1}".format(_InputData.max_peak_cut, bpm_with_spike_counter)
        print "Startturn: {0} Maxturns: {1}".format(_InputData.startturn_human, _InputData.maxturns_human)
        print "Number of turns: " + str(self.number_of_turns)
        
        print "Horizontal BPMs: " + str(self.dictionary_plane_to_bpms[PLANE_X].get_number_of_bpms())
        print "Vertical   BPMs: "   + str(self.dictionary_plane_to_bpms[PLANE_Y].get_number_of_bpms())
        
        bpms_removed = bpm_with_exact_zero_counter + known_bad_bpms_counter + flatbpm_counter + bpm_with_spike_counter
        good_bpms = self.dictionary_plane_to_bpms[PLANE_X].get_number_of_bpms() + self.dictionary_plane_to_bpms[PLANE_Y].get_number_of_bpms()
        total_bpms = bpms_removed + good_bpms
        print "(Statistics for file reading) Total BPMs: {0}, Good BPMs: {1} ({2:2.2f}%), bad BPMs: {3} ({4:2.2f}%)"\
            .format(total_bpms, good_bpms, 100.0*good_bpms/total_bpms, bpms_removed, 100.0*bpms_removed/total_bpms)
        
        print "Reading done"
        if PRINT_TIMES:
            print ">>Time for init (read file): {0}s".format(time.time() - time_start)

    def remove_bpms(self, badbpm_indices, plane, reason):
        """Removes the given BPMs for the given plane with the given reason"""
        badbpms = self.dictionary_plane_to_bpms[plane].remove_badbpm_and_get_badbpms(badbpm_indices, reason)
        self.bad_bpm_file.add_badbpms(badbpms)

    def write_file(self):
        """Writes the new file"""
        time_start = time.time()
        self.bad_bpm_file.header = self.header
        self.bad_bpm_file.write_badbpms()
        temp_path = self.path_to_sddsfile + ".tmp_svd_clean"
        print "writing data to temp file: " + temp_path

        with open(temp_path, "w") as temp_file:
            temp_file.write(self.header)
            if "NTURNS calculated" not in self.header:
                temp_file.write("#NTURNS calculated: " + str(self.number_of_turns) + "\n")
            if not _InputData.use_test_mode:
                temp_file.write("#Modified: {0} By: {1}\n".format(time.strftime("%Y-%m-%d#%H:%M:%S"), __file__))  # to get only the filename: os.path.basename(__file__)
            self.write_bpm_data_to_file(PLANE_X, temp_file)
            self.write_bpm_data_to_file(PLANE_Y, temp_file)

        print "writing done. rename to: " + str(_InputData.outputfile)
        if not _InputData.use_test_mode:
            try:
                os.remove(_InputData.outputfile)
            except OSError:
                pass  # Will be raised if outputfile does not exist or no access -vimaier
            os.rename(temp_path, _InputData.outputfile)
        
        if PRINT_TIMES:
            print ">> Time for write_file: {0}s".format(time.time() - time_start)

    def write_bpm_data_to_file(self, plane, file_to_write):
        """Writes the BPM data for the given plane"""
        current_bpms = self.dictionary_plane_to_bpms[plane]
        for i, bpm_name_location_plane_item in enumerate(current_bpms.bpm_name_location_plane):
            file_to_write.write("{0} {1[0]:<15} {1[1]:>12.5f} ".format(plane, bpm_name_location_plane_item))
            # {0} is the first argument, {1[0]} from the second argument the first entry, {1[1] from the second argument the second entry -tbach
            # :>15 means right aligned, filled up to 15 characters. example: "   BPMYA.4L1.B1" -tbach
            # :>12.5 means right aligned float, filled up to 12 characters and fixed precision of 5. Example: " 23347.14262" -tbach
            current_bpms.bpm_data[i].tofile(file_to_write, sep=" ", format="% 8.5f")
            # % 8.5f is a float, filled up to 8 characters and fixed precision of 5. if negative, preceded by a sign, if positive by a space -tbach
            # Example1: " -0.94424", example2: "  1.25630" -tbach
            file_to_write.write("\n")


class _BadBpmFile(object):
    """ Represents a bad bpm file
    Handles interaction with the bad bpm file """

    def __init__(self, path_to_sddsfile):
        self.path_to_sddsfile = path_to_sddsfile
        self.path_to_badbpmfile = path_to_sddsfile.replace(".new", "") + ".bad"
        self.header = ""
        self.lines_to_write = []

    def __write_header(self, filehandle_badbpm):
        """ Writes the header for the bad bpm file """
        filehandle_badbpm.write(self.header)
        if not _InputData.use_test_mode:
            filehandle_badbpm.write("#Modified: {0} By: {1}\n".format(time.strftime("%Y-%m-%d#%H:%M:%S"), __file__))  # to get only the filename: os.path.basename(__file__)
        filehandle_badbpm.write("#Bad BPMs from: {0}\n".format(self.path_to_sddsfile))

    def add_badbpm(self, badbpm):
        """ Adds a new bad BPM to write """
        self.lines_to_write.append(badbpm.to_string)

    def add_badbpms(self, list_badbpm):
        """Adds all bad BPM to write """
        for badbpm in list_badbpm:
            self.add_badbpm(badbpm)

    def write_badbpms(self):
        """Writes all the current bad BPM"""
        if not self.lines_to_write:
            return
        time_start = time.time()
        print "Create file for bad BPMs: ", self.path_to_badbpmfile

        with open(self.path_to_badbpmfile, "w") as file_badbpms:
            self.__write_header(file_badbpms)
            for line in self.lines_to_write:
                file_badbpms.write(line)
        
        if PRINT_TIMES:
            print ">> Time for write_badbpms: {0}s".format(time.time() - time_start)


class _Bpms(object):
    """This represents some BPMs. It is used to distinguish between BPMs for X and Y"""

    def __init__(self):
        self.bpm_data = []
        self.bpm_name_location_plane = []

    def get_number_of_bpms(self):
        """ Get the number of BPMs """
        return len(self.bpm_name_location_plane)

    def remove_badbpm_and_get_badbpms(self, badbpm_indices, reason):
        """Removes the given BPM indices from the current BPM list and returns a list of the given bad BPMs"""
        badbpms = []
        for index in sorted(badbpm_indices, reverse=True):
            badbpms.append(_BadBpm(self.bpm_name_location_plane.pop(index), self.bpm_data.pop(index), reason))
        return badbpms


class _BadBpm(object):
    """Represents a bad BPM"""
    
    def __init__(self, *args, **kwargs): #this is python, we can not have multiple constructors :/ -tbach
        if len(args) == 3:
            self.init_with_variables(args[0], args[1], args[2])
        elif len(args) == 2:
            self.init_with_line(args[0], args[1])
        else:
            print "WARNING. unknown arguments, check source. len(args): " + str(len(args))
            self.init_with_line(str(args), "unknown arguments")

    def init_with_variables(self, bpm_name_location_plane, data, reason):
        #print "badbpm " + str(reason) #show all removed bad bpms with reason -tbach
        name = bpm_name_location_plane[0]
        location = bpm_name_location_plane[1]
        plane = bpm_name_location_plane[2]
        
        self.line_to_write = "{0} {1:<15} {2:>12.5f} {3}\n#{4}\n" \
            .format(plane, name, location, " ".join(["{0:>12.5f}".format(x) for x in data]), reason)
        # for string formatting explanation, see write_bpm_data_to_file -tbach
        
    def init_with_line(self, line, reason):
        self.line_to_write = "{0}#{1}\n".format(line, reason)

    @property
    def to_string(self):
        """ returns the String representation for this bad bpm"""
        return self.line_to_write


class _SvdHandler(object):
    """ Main workload class. Handles the SVD"""

    def __init__(self):
        self.sddsfile = _SddsFile(_InputData.inputfile)
        self.sddsfile.init()
        print ""

        time_start = time.time()
        self.do_svd_clean(PLANE_X)
        self.do_svd_clean(PLANE_Y)
        if PRINT_TIMES:
            print ">> Time for svdClean (X&Y): {0}s".format(time.time() - time_start)
        print ""

        self.sddsfile.write_file()
        print ""

    def do_svd_clean(self, plane):
        """Does a SVD clean on the given plane"""
        time_start = time.time()
        print "(plane {0}) removing noise floor with SVD".format(plane)

        A = numpy.array(self.sddsfile.dictionary_plane_to_bpms[plane].bpm_data)
        number_of_bpms = A.shape[0]

        if number_of_bpms <= 10:
            sys.exit("Number of bpms <= 10")

        # normalise the matrix
        sqrt_number_of_turns = numpy.sqrt(A.shape[1])
        A_mean = numpy.mean(A)
        A = (A - A_mean) / sqrt_number_of_turns

        if PRINT_TIMES:
            print ">> Time for svdClean (before SVD call): {0}s".format(time.time() - time_start)
        USV = self.get_singular_value_decomposition(A)
        if PRINT_TIMES:
            print ">> Time for svdClean (after SVD call): {0}s".format(time.time() - time_start)

        # remove bad BPM by SVD -tbach
        good_bpm_indices = self.remove_bad_bpms_and_get_good_bpm_indices(USV, plane)
        USV = (USV[0][good_bpm_indices], USV[1], USV[2])
        number_of_bpms = len(good_bpm_indices)

        #----SVD cut for noise floor
        if _InputData.singular_values_amount_to_keep < number_of_bpms:
            print "(plane {0}) amount of singular values to keep: {1}".format(plane, _InputData.singular_values_amount_to_keep)
            USV[1][_InputData.singular_values_amount_to_keep:] = 0
        else:
            print "requested more singular values than available(={0})".format(number_of_bpms)

        A = matrixmultiply(USV[0], matrixmultiply(numpy.diag(USV[1]), USV[2]))
        # A0 * (A1 * A2) should require less operations than (A0 * A1) * A2,
        #  because of the different sizes
        # A0 has (M, K), A1 has (K, K) and A2 has (K, N) with K=min(M,N)
        # Most of the time, the number of turns is greater then
        #  the number of BPM, so M > N
        # --tbach
        A = (A * sqrt_number_of_turns) + A_mean
        bpmres = numpy.mean(numpy.std(A - self.sddsfile.dictionary_plane_to_bpms[plane].bpm_data, axis=1))
        print "(plane {0}) Average BPM resolution: ".format(plane) + str(bpmres)
        self.sddsfile.dictionary_plane_to_bpms[plane].bpm_data = A
        if PRINT_TIMES:
            print ">> Time for do_svd_clean: {0}s".format(time.time() - time_start)
            

    def get_singular_value_decomposition(self, matrix):
        """Calls the SVD, returns a USV representation
        For details, see numpy documentation"""
        return numpy.linalg.svd(matrix, full_matrices=False)  # full matrices do not have any interesting value for us -tbach

    def remove_bad_bpms_and_get_good_bpm_indices(self, USV, plane):
        """ Removes bad(dominant) BPMs and gets the indices from good BPMs for the given plane"""
        time_start = time.time()

        # A is [u,s,v] with u * np.diag(s) * v = original matrix -tbach
        U_t_abs = numpy.transpose(abs(USV[0]))  # This creates a view, which is nice and fast -tbach
        # What happens here?
        # From the SDDS ASCII file, We have a BPM x Turns matrix.
        # Let B = Number of BPM, T = Number of Turns, then matrix size is (B,T)
        # If we do the SVD, we have U,S,V. U has the size (B,x) and V (x,T)
        # If we transpose, U^t has the size (x,B)
        # Now, we can look in every row to get the maximum value for every BPM.
        # If one BPM is dominating, we remove it as a bad BPM
        # --tbach

        badbpm_indices = set()  # we do not want duplicates -tbach

        for row_index in range(len(U_t_abs)):
            max_index = numpy.argmax(U_t_abs[row_index])
            max_value = U_t_abs[row_index][max_index]
            if (max_value > _InputData.single_svd_bpm_threshold):
                badbpm_indices.add(max_index)

        if PRINT_DEBUG:
            print "Bad BPM indices: " + str(badbpm_indices)

        number_of_badbpms = len(badbpm_indices)
        if number_of_badbpms > 0:
            print "(plane {0}) Bad BPMs from SVD detected. Number of BPMs removed: {1}".format(plane, len(badbpm_indices))

        # add the bad BPMs to the general list of bad BPMs -tbach
        reason_for_badbpm = "Detected from SVD, single peak value is greater then " + str(_InputData.single_svd_bpm_threshold)
        self.sddsfile.remove_bpms(badbpm_indices, plane, reason_for_badbpm)

        number_of_bpms = USV[0].shape[0]
        goodbpm_indices = range(number_of_bpms)
        for value in badbpm_indices:
            goodbpm_indices.remove(value)

        if PRINT_TIMES:
            print ">> Time for removeBadBpms: {0}s".format(time.time() - time_start)
        return goodbpm_indices


#===================================================================================================
# main invocation
#===================================================================================================
def _start():
    """Starter function to avoid polluting global namespace with variable 'options'."""
    options = _parse_args()
    clean_sdds_file(options)

if __name__ == "__main__":
    _start()
