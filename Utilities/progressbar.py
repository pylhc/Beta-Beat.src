'''
Created on Aug 9, 2016

@author: awegsche
'''

import sys
from time import time

PROGRESSBAR_LENGTH = 40


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
    sys.stdout.write(" [\33[32m" + "#" * pos + "\33[0m" + "#" * (PROGRESSBAR_LENGTH - pos) + "] " + "{0:5.1f} %".format(x) + chr(8) * (PROGRESSBAR_LENGTH + 11))
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
