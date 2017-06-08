#
# CODE TO ANALYSE AN IR NONLINEARITY SCAN
#
import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD/")
import string
import os
import re
import math
import csv
#from numpy import *
import numpy as np
from optparse import OptionParser
import atexit
from operator import itemgetter
import datetime as dt
import matplotlib
#matplotlib.use('Qt4gg')
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
import scipy.odr
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.optimize import leastsq
from scipy.optimize import basinhopping
import multiprocessing as mp
import pytimber
import shutil
#import progressbar
import time as ti
import copy as cp

def parse_timestamp(thistime):
    accepted_time_input_format = ['%Y-%m-%d %H:%M:%S.%f','%Y-%m-%d %H:%M:%S','%Y-%m-%d_%H:%M:%S.%f','%Y-%m-%d_%H:%M:%S']
    for fmat in accepted_time_input_format:
        try:
            dtobject=dt.datetime.strptime(thistime,fmat)
            return dtobject
        except ValueError:
            pass
    timefmatstring=''
    for fmat in accepted_time_input_format:
        timefmatstring=timefmatstring+'\"'+fmat+'\" ,   '
    sys.tracebacklimit = 0
    raise ValueError('No appropriate input format found for start time of scan (-s).\n ---> Accepted input formats are:   '+timefmatstring)
###########################################################################################################################################################################
def convert_name(name):
    new_name = re.sub("\s+","_", name)
    new_name = re.sub("/","_",name)
    return new_name
#### FUNCTION TO CONVERT A DATETIME OBJECT TO THE TIMESTAMP FORMAT REQUIRED BY THE ONLINE MODEL. SHOULD CLARIFY UTC POINT HERE
def convert_to_OM_input_format(dtobject):
    OM_timestamp = dtobject.strftime('%Y-%m-%d %H:%M:%S.%f')
    return OM_timestamp
def convert_to_timber_input_format(dtobject): ### for now this is identical to OM input format. but separated the functions in case of future differences.
    timber_timestamp = dtobject.strftime('%Y-%m-%d %H:%M:%S.%f')
    return timber_timestamp
def convert_to_data_output_format(dtobject): ### function to help write output from datetime objects in standard format throughout code
    output_timestamp=dtobject.strftime('%Y-%m-%d %H:%M:%S.%f')
    return output_timestamp
def convert_from_data_output_format(output_timestamp): ### function to help read in times from existing output, generated with the function above.
    dtobject=dt.datetime.strptime(output_timestamp,'%Y-%m-%d %H:%M:%S.%f')
    return dtobject
###########################################################################################################################################################################
def unix_to_dtime_converter(unixtimestamp):
    epoch='1970-01-01 00:00:00.000'
    epochtime=dt.datetime.strptime(epoch,'%Y-%m-%d %H:%M:%S.%f')
    unixtime=float(unixtimestamp) 
    dtime=dt.timedelta(0,unixtime,0)
    dtime=epochtime+dtime
    return dtime
def dtime_to_unix_converter(dtobject):
    epoch='1970-01-01 00:00:00.000'
    epochtime=dt.datetime.strptime(epoch,'%Y-%m-%d %H:%M:%S.%f')
    timediff=(dtobject-epochtime)
    nsec=(timediff.microseconds / 1e6 + timediff.seconds + timediff.days*24*60*60)
    unixtimestamp=nsec
    return unixtimestamp
def local_to_utc_converter(localtime):
    utctime=''
    return utctime
def utc_to_local_converter(utctime):
    localtime=''
    return localtime



######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################




def parse_options():
    parser = OptionParser()
    ####################################################################
    parser.add_option("-k", "--knob1",
                      help="Name of the knob which has been scanned for Beam1. By default this is the knob for both beams.",
                      metavar="KNOB1", default="LHCBEAM/IP5-XING-H-MURAD", dest="knob1")
    parser.add_option("-j", "--knob2",
                      help="Name of the knob which has been scanned for Beam2. Use if have knobs like LHCBEAM1/ and LHCBEAM2/ rather thank LHCBEAM/. If have LHCBEAM/ type knob, only specify option \'-k\'.",
                      metavar="KNOB2", default="LHCBEAM/IP5-XING-H-MURAD", dest="knob2")
    ####################################################################
    parser.add_option("-o", "--outputpath",
                      help="Path to which results and extracted data are output",
                      metavar="OUTPUTPATH", default="./", dest="outputpath")
    parser.add_option("-f", "--datafile",
                      help="Path to file with platteau times",
                      metavar="DATAFILE", default="./accepted_platteaus.dat", dest="datafile")
    ###################################################################
    (options, args) = parser.parse_args()
    knob1=options.knob1
    knob2=options.knob2
    outputpath=options.outputpath
    outputpath=os.path.abspath(outputpath)+'/'
    datafile=options.datafile
    datafile=os.path.abspath(datafile)
    ###################################################################
    ###### SANITY CHECKS OF SPECIFIED OPTIONS #########################
    #############
    if not os.path.exists(outputpath):
        sys.exit('Could not find specified output path: '+outputpath)
    if not os.path.exists(datafile):
            sys.exit('Could not find specified input file: '+datafile)
    ###################################################################
    return knob1,knob2,outputpath,datafile



def getknobsetting(plattime,knob,outputpath):
    timestamp_OM = convert_to_OM_input_format(plattime)
    print '---> EXTRACTING KNOB SETTINGS FOR KNOB:'+knob+' at '+timestamp_OM

    command='/afs/cern.ch/eng/lhc_online_model/dev/bin/lhc-model-extractor.dev_knobs.sh -beam B1 -time \"'+timestamp_OM+'\" -Knames \"'+knob+'\" -Ktype KNOBVALUE -p '+outputpath
    #os.system(command)
    os.popen(command)
    OMoutputfile='knobs.madx'
    OMoutput=open(outputpath+OMoutputfile,'r')
    for line in OMoutput.readlines():
        if re.search('beamprocess',line):
            beam_process = line.partition('beamprocess')[2].partition('at ')[0].strip()
            print '---> Beam Process = '+beam_process
        elif re.search('! Extract Delta values for knob',line):
            filetime = line.partition('at time')[2].strip().strip('.').strip()
            dtfiletime = parse_timestamp(filetime)
            print '---> Last trim of knob was at: '+filetime
        elif re.search('! Value of knob',line):
            knob_setting = float(line.partition('=')[2].partition(';')[0].strip())
            knob_trim_var = line.partition('=')[0].strip() # get the variable name for the knob definition. will be used to check for circuits defined at zero in next step.
            print '---> Knob setting = '+str(knob_setting)
            print '\n'
    OMoutput.close()
    return knob_setting,dtfiletime
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################

def readplatteaus(datafile,knob1,knob2,outputpath):
    print 'Extracting Trim History for platteaus defined in: '+datafile+'\n\n'
    rplat=open(datafile,'r')
    csvrplat=csv.reader(rplat,delimiter=',',skipinitialspace=True)

    knobtrimvalues=[]
    for row in csvrplat:
        if re.search('B1_min',row[1]):
            continue
        else:
            platteau_index=int(row[0])
            b1_plat_start=row[1].replace('+00:00','')
            b1_plat_end=row[2].replace('+00:00','')
            b2_plat_start=row[3].replace('+00:00','')
            b2_plat_end=row[3].replace('+00:00','')
            dt_b1start=parse_timestamp(b1_plat_start)
            dt_b1end=parse_timestamp(b1_plat_end)
            dt_b2start=parse_timestamp(b2_plat_start)
            dt_b2end=parse_timestamp(b2_plat_end)

            print '\n\n\n'
            print 'EXTRACTING KNOB SETTING FROM PLATTEAU: '+str(platteau_index)
            print '------------------------------------------------------------'
            print 'BEAM1 EXTRACTION:'
            print '-----------------\n'
            b1start_knobsetting,dt_OMb1starttime = getknobsetting(dt_b1start,knob1,outputpath)
            b1end_knobsetting,dt_OMb1endtime     = getknobsetting(dt_b1end,knob1,outputpath)
            if b1start_knobsetting!=b1end_knobsetting:
                print ''
                endtimediff=dt_b1end-dt_OMb1endtime
                endtimediff=endtimediff.total_seconds()
                print 'WARNING WARNING WARNING:    KNOB VALUES AT START AND END OF BEAM1 PLATTEAU ARE NOT THE SAME'
                print 'WARNING WARNING WARNING:    TIME BETWEEN PLATTEAUFILE END, AND KNOB TRIM IS '+str(endtimediff)+'s'
                print 'WARNING WARNING WARNING:    IF SMALL MAYBE JUST AN ISSUE WITH TIMBER VS LSA LOGGING TIMES'
                print 'WARNING WARNING WARNING:    CONTINUING ON ASSUMPTION YOU KNOW WHAT YOU ARE DOING!?! VALUES TAKEN ARE FROM THE STARTTIME OF THE PLATTEAU'
                print ''

            print '\n'
            print 'BEAM2 EXTRACTION:'
            print '-----------------\n'
            b2start_knobsetting,dt_OMb2starttime = getknobsetting(dt_b2start,knob2,outputpath)
            b2end_knobsetting,dt_OMb2endtime     = getknobsetting(dt_b2end,knob2,outputpath)
            if b2start_knobsetting!=b2end_knobsetting:
                print ''
                endtimediff=dt_b2end-dt_OMb2endtime
                endtimediff=endtimediff.total_seconds()
                print 'WARNING WARNING WARNING:    KNOB VALUES AT START AND END OF BEAM2 PLATTEAU ARE NOT THE SAME'
                print 'WARNING WARNING WARNING:    TIME BETWEEN PLATTEAUFILE END, AND KNOB TRIM IS '+str(endtimediff)+'s'
                print 'WARNING WARNING WARNING:    IF SMALL MAYBE JUST AN ISSUE WITH TIMBER VS LSA LOGGING TIMES'
                print 'WARNING WARNING WARNING:    CONTINUING ON ASSUMPTION YOU KNOW WHAT YOU ARE DOING!?! VALUES TAKEN ARE FROM THE STARTTIME OF THE PLATTEAU'
                print ''

            knobdata=[platteau_index,b1start_knobsetting,b2start_knobsetting]
            knobtrimvalues.append(knobdata)
    rplat.close()

    return knobtrimvalues

######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################





if __name__ == "__main__":
    ####################################################################
    print '\n\n'
    print '----------------------------------------------------------------------------------'
    print 'A code to obtain and analyse IR-orbit scan data'
    print '----------------------------------------------------------------------------------'
    print 'Need annaconda to run: /afs/cern.ch/work/o/omc/anaconda/bin/python'
    print 'Running version: '+str(sys.version)
    print '\n'
    ####################################################################
    ## PROVIDE GENERIC VALUES FOR THE CODE, EG UTC/LOCAL CONVERSION, EPOCH DATETIME OBJECT
    utcconvert=dt.timedelta(0,7200)
    ####################################################################
    ## PARSE THE OPTIONS FOR THE CODE
    knob1,knob2,outputpath,datafile = parse_options()

    knobtrimvalues = readplatteaus(datafile,knob1,knob2,outputpath)

    wout=open(outputpath+'data.platteau_knob_settings.csv','w')
    csvwout=csv.writer(wout,delimiter=',')
    for p in range(len(knobtrimvalues)):
        csvwout.writerow(knobtrimvalues[p])
    wout.close()
