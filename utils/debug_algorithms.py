# uncompyle6 version 3.2.3
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.5.4 |Anaconda 4.2.0 (64-bit)| (default, Nov 20 2017, 18:44:38) 
# [GCC 7.2.0]
# Embedded file name: /media/awegsche/HDD/CernBox/Beta-Beat.src/GetLLM/algorithms/debug_algorithms.py
# Compiled at: 2018-02-14 12:23:46
from struct import pack
from utils import logging_tools
LOGGER = logging_tools.get_logger(__name__)
DOUBLE_T = 3
BPM_T = 6
MATRIX_T = 7
COMBINATION_T = 5
# TODO: this should eventually be put in some kind of class. But for now, meh.
DFILE = None

def create_debugfile(filename):
    global DFILE
    if DFILE is not None:
        raise IOError("The debugfile has already been created. Please save and cloase before "
                      "creating a new one.")
    LOGGER.debug(("creating debugfile '{}'").format(filename))
    DFILE = open(filename, 'wb')
    DFILE.write('BDBG')


def close_file():
    global DFILE
    LOGGER.debug('Closing debugfile')
    DFILE.write('"')
    DFILE.close()
    DFILE = None


def start_write_bpm(bpm_name, s, beta, alpha, mu):
    global DFILE
    DFILE.write(pack(('=bi{}sdddd').format(len(bpm_name)), BPM_T, len(bpm_name),
                     bpm_name, s, beta, alpha, mu))


def start_write_combinations(count):
    global DFILE
    DFILE.write(pack('=bi', COMBINATION_T, count))


def write_end():
    global DFILE
    DFILE.write('!')


def write_bpm_combination(n1, n2, beta, weight):
    global DFILE
    DFILE.write(pack('=iidd', n1, n2, beta, weight))


def write_matrix(matrix, name):
    global DFILE
    n, m = matrix.shape
    DFILE.write(pack(('=bi{}sii').format(len(name)), MATRIX_T, len(name), name, n, m))
    for y in range(m):
        for x in range(n):
            DFILE.write(pack('d', matrix[(x, y)]))


def write_double(name, value):
    global DFILE
    DFILE.write(pack(('=bi{}sd').format(len(name)), DOUBLE_T, len(name), name, value))
