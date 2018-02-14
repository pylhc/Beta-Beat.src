import sys
import os
from struct import pack
from utils import logging_tools

LOGGER = logging_tools.get_logger(__name__)

DOUBLE_T = 3
BPM_T = 6
MATRIX_T = 7
COMBINATION_T = 5

dfile = None


def create_debugfile(filename):
    global dfile
    LOGGER.debug("creating debugfile '{}'".format(filename))
    dfile = open(filename, "wb")
    dfile.write("BDBG")

def close_file():
    global dfile
    LOGGER.debug("CLosing debugfile")
    dfile.write("\x22")
    dfile.close()

def start_write_bpm(bpm_name, s, beta, alpha, mu):
    global dfile
    dfile.write(pack("=bi{}sdddd".format(len(bpm_name)), BPM_T, len(bpm_name), bpm_name, s, beta, alpha, mu))

def start_write_combinations(count):
    global dfile
    dfile.write(pack("=bi", COMBINATION_T, count))

def write_end():
    global dfile
    dfile.write("\x21")

def write_bpm_combination(n1, n2, beta, weight):
    global dfile
    dfile.write(pack("=iidd", n1, n2, beta, weight))

def write_matrix(matrix, name):
    global dfile
    (n,m) = matrix.shape

    dfile.write(pack("=bi{}sii".format(len(name)), MATRIX_T, len(name), name, n, m))

    for y in range(m):
        for x in range(n):
            dfile.write(pack("d", matrix[x,y]))

def write_double(name, value):
    global dfile
    dfile.write(pack("=bi{}sd".format(len(name)), DOUBLE_T, len(name), name, value))
