'''
Created on Aug 9, 2016

@author: awegsche
'''

import sys
from time import time

from numpy.random import random

PROGRESSBAR_LENGTH = 40
PROTON_PATHLENGTH = 17

def startProgress(title):
    global progress_x
    sys.stdout.write(" [" + "#" * PROGRESSBAR_LENGTH + "]" + chr(8) * (PROGRESSBAR_LENGTH + 3))
    sys.stdout.flush()
    progress_x = 0


def startProgress_(title):
    sys.stdout.write(" [" + "#" * PROGRESSBAR_LENGTH + "]" + chr(8) * (PROGRESSBAR_LENGTH + 1))
    sys.stdout.flush()


def progress(x):
    pos = int(x * PROGRESSBAR_LENGTH)
    sys.stdout.write(" [\33[32m" + "#" * pos + "\33[0m" + "#" * (PROGRESSBAR_LENGTH - pos) + "] " + "{0:5.1f} %".format(x*10) + chr(8) * (PROGRESSBAR_LENGTH + 11))
    sys.stdout.flush()
    

def endProgress():
    sys.stdout.write(" [\33[32m" + "#" * PROGRESSBAR_LENGTH + "\33[0m] 100.0 %\n")
    sys.stdout.flush()
    
    
def startProgressBW(title):
    global progress_x
    sys.stdout.write(title + ": [" + " " * PROGRESSBAR_LENGTH + "]" + chr(8) * (PROGRESSBAR_LENGTH + 1))
    sys.stdout.flush()
    progress_x = 0


def progressBW(x):
    global progress_x
    x = int(x * PROGRESSBAR_LENGTH // 100)
    sys.stdout.write("#" * (x - progress_x))
    sys.stdout.flush()
    progress_x = x


def endProgressBW():
    sys.stdout.write("#" * (PROGRESSBAR_LENGTH - progress_x))
    sys.stdout.flush()


def startTime():
    global __time_started__
    global __step__
    __time_started__ = time()
    __step__ = 0
    
    
def printStep(step):
    global __time_started__
    position = step % (PROGRESSBAR_LENGTH - 3)
    
    elapsed = time() - __time_started__
    elapsed_str = ""
    size = 0
    
    if elapsed < 60:
        elapsed_str = "{0:6.1f}s".format(elapsed)
        size = 7
    else:
        min_part = int(elapsed // 60)
        sec_part = elapsed - (float)(min_part * 60.0)
        elapsed_str = "{0:2d}:{1:04.1f}".format(min_part, sec_part)
        size = 7
    
    sys.stdout.write(" [" + "=" * (position) + "\33[32m===\33[0m" + "=" * (PROGRESSBAR_LENGTH - position - 3) + "] " +
                     elapsed_str + chr(8) * (PROGRESSBAR_LENGTH + 5 + size))
    sys.stdout.flush()
   
    
def printStepPercent(percent):
    global __time_started__
    global __step__
    position = __step__ % (PROGRESSBAR_LENGTH - 3)
    __step__ += 1
    
    elapsed = time() - __time_started__
    elapsed_str = ""
    size = 0
    
    if elapsed < 60:
        elapsed_str = "{0:6.1f}s".format(elapsed)
        size = 7
    else:
        min_part = int(elapsed // 60)
        sec_part = elapsed - (float)(min_part * 60.0)
        elapsed_str = "{0:2d}:{1:04.1f}".format(min_part, sec_part)
        size = 7
    
    sys.stdout.write(" [" + "=" * (position) + "\33[32m===\33[0m" + "=" * (PROGRESSBAR_LENGTH - position - 3) + "] " +
                     elapsed_str + " {0:4.1f} %".format(percent) + chr(8) * (PROGRESSBAR_LENGTH + 12 + size))
    sys.stdout.flush()
    
    
def stopPrint():
    global __time_started__
    elapsed = time() - __time_started__
    elapsed_str = ""
    
    if elapsed < 60:
        elapsed_str = "{0:6.1f}s".format(elapsed)
    else:
        min_part = int(elapsed // 60)
        sec_part = elapsed - (float)(min_part * 60.0)
        elapsed_str = "{0:2d}:{1:04.1f}".format(min_part, sec_part)

    sys.stdout.write(" [" + "=" * PROGRESSBAR_LENGTH + "] " +
                     elapsed_str + "\n")
    sys.stdout.flush()


def setupProtons():
    global  __step__
    __step__ = 0


def progressProtons():
    print_protons()
    sys.stdout.flush()


def print_protons():
    global __step__
    position = __step__ % (PROTON_PATHLENGTH + 5)
    if position < PROTON_PATHLENGTH:
        sys.stdout.write("[" + " " * (position) + "\33[38;2;90;90;255mp+\33[0m" +
                         " " * (2 * (PROTON_PATHLENGTH - position) + 1) + "\33[38;2;90;90;255mp+\33[0m" + " " * (
                         position) + "] " +
                         chr(8) * (2 * PROTON_PATHLENGTH + 8))
    elif position % 2 == 1:
        sys.stdout.write(
            "[" + " " * (PROTON_PATHLENGTH + 2) + "\33[38;2;200;200;50mX\33[0m" + " " * (PROTON_PATHLENGTH + 2) + "]" +
            chr(8) * (2 * PROTON_PATHLENGTH + 8))
    else:
        sys.stdout.write(
            "[" + " " * (PROTON_PATHLENGTH + 2) + "\33[38;2;200;0;50mX\33[0m" + " " * (PROTON_PATHLENGTH + 2) + "]" +
            chr(8) * (2 * PROTON_PATHLENGTH + 8))
    __step__ = __step__ + 1


def setupLumi(seed=None):
    setupProtons()
    global __iLumi__
    __iLumi__ = random()*150+100

def progressLumi(percent):
    text = "[Integrated luminosity: {0:5.1f} / {1:5.1f}]                  \n".format(percent * __iLumi__, __iLumi__)
    sys.stdout.write(text)
    print_protons()
    sys.stdout.write(chr(8) * len(text))
    sys.stdout.flush()

def endLumi():
    text = "[Integrated luminosity: {0:5.1f} / {0:5.1f}]                  \n".format(__iLumi__)
    sys.stdout.write(text)
    print_protons()
    sys.stdout.write(chr(8) * len(text))
    sys.stdout.flush()