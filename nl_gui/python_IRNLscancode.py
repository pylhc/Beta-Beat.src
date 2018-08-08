###########################################################
##### CODE TO ANALYSE AN IR NONLINEARITY SCAN
###########################################################
#
### EHM 28/04 editing headers to match felix's gui input
#
#
#
import sys
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

######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
###########################################################################################################################################################################
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

###########################################################################################################################################################################
#### READ KNOBNAME, START AND END TIME OF SCAN, FLAG TO EXTRACT DATA / PATH TO PREVIOUSLY EXTRACTED DATA, FLAG TO CHECK AC DIPOLE KICKS
def parse_options():
    parser = OptionParser()
    ####################################################################
    parser.add_option("-k", "--knob1",
                      help="Name of the knob which has been scanned for Beam1. By default this is the knob for both beams.",
                      metavar="KNOB1", default="LHCBEAM/IP5-XING-H-MURAD", dest="knob1")
    parser.add_option("-j", "--knob2",
                      help="Name of the knob which has been scanned for Beam2. Use if have knobs like LHCBEAM1/ and LHCBEAM2/ rather thank LHCBEAM/. If have LHCBEAM/ type knob, only specify option \'-k\'.",
                      metavar="KNOB2", default=None, dest="knob2")
    ####################################################################
    parser.add_option("-s", "--starttime",
                      help="Start time of the scan in local time",
                      metavar="STARTTIME", default="2016-08-22 10:45:00.000", dest="starttime")
    parser.add_option("-e", "--endtime",
                      help="End time of the scan in local time",
                      metavar="ENDTIME", default="2016-08-22 11:38:00.000", dest="endtime")
    parser.add_option("-u", "--utcinput",action="store_true",
                      help='Flag for UTC time input. If present input times for the -s and -e option are assumed to be in UTC, if absent Geneva local time. Not that default values for start and end time test case are local.',
                      metavar="UTCINPUT", default=False, dest="utcinput")
    ####################################################################
    parser.add_option("-c", "--cofb", action="store_true",
                  help="Flag to specify whether an automatic closed orbit feed-back is present. If true will NOT perform checks to exclude individual MCB trims within the knob platteaus.",
                  metavar="OFB", default=False, dest="ofb")
    parser.add_option("-a", "--acdipole", action="store_true",
                  help="Flag to specify if AC-dipole kicks were performed during the scan",
                  metavar="ACDIPOLE", default=False, dest="acdipole")
    ###################################################################
    parser.add_option("-x", "--extract", action="store_true",
                      help="Flag for whether to extract data from timber (-x), or look for previously extracted data",
                      metavar="EXTRACT", default=False, dest="extract")
    parser.add_option("-o", "--outputpath",
                      help="Path to which results and extracted data are output",
                      metavar="OUTPUTPATH", default="./", dest="outputpath")
    parser.add_option("-d", "--datapath",
                      help="Path to previously extracted data, expected to be None if -x is true",
                      metavar="DATAPATH", default=None, dest="datapath")
    ###################################################################
    (options, args) = parser.parse_args()
    knob1=options.knob1
    knob2=options.knob2
    knobs=[]
    knobs.append(knob1)
    if knob2!=None:
        knobs.append(knob2)
    starttime=options.starttime
    endtime=options.endtime
    utcinput=options.utcinput
    OFB=options.ofb
    acdipole=options.acdipole
    extract=options.extract
    outputpath=options.outputpath
    datapath=options.datapath
    ###################################################################
    dtstart = parse_timestamp(starttime)
    dtend = parse_timestamp(endtime)
    ###################################################################
    outputpath=os.path.abspath(outputpath)+'/'
    if datapath!=None:
        datapath=os.path.abspath(datapath)+'/'
    ###################################################################
    ###### SANITY CHECKS OF SPECIFIED OPTIONS #########################
    #############
    if knob1==knob2:
        sys.exit('Specified same knob for \'-k\' and \'-j\' options. If have LHCBEAM/ type knob acting on both, specify only \'-k\'.')
    if len(knobs)>2:
        sys.exit('Why are there too many knobs!?!')
    for i in range(len(knobs)):
        if re.search('^LHCBEAM/',knobs[i]) or re.search('^LHCBEAM1/',knobs[i]) or re.search('^LHCBEAM2/',knobs[i]):
            pass
        else:
            sys.exit('Knob name provided: ['+knobs[i]+']...\n ...do not include the LHCBEAM specifier. This is required to distinguish 2-beam knobs, and beam1/2 only knobs.')
    if len(knobs)==1:
        pass
    elif len(knobs)>1:
        if re.search('^LHCBEAM1/',knobs[0]) and re.search('^LHCBEAM1/',knobs[1]):
            sys.exit('Provided 2 knob names: '+knobs[0]+'\n                       '+knobs[1]+'...\n ...but both specify LHCBEAM1/')
        elif re.search('^LHCBEAM2/',knobs[0]) and re.search('^LHCBEAM2/',knobs[1]):
            sys.exit('Provided 2 knob names: '+knobs[0]+'\n                       '+knobs[1]+'...\n ...but both specify LHCBEAM2/')
    #############
    if dtstart>=dtend:
        sys.exit('Start time ('+starttime+') must be before end time ('+endtime+')')
    #############
    if not os.path.exists(outputpath):
        sys.exit('Could not find specified output path: '+outputpath)
    #############
    if extract==True:
        if datapath!=None:
            sys.exit('Have specified to extract data from timber (option \'-x\'), but also provided a path to existing data (option \'-d\')???')
    elif  extract==False:
        if datapath==None:
            sys.exit('Extraction from timber (\'-X\') not requested, but no path to existing data (\'-d\') provided. Try again.')
        if not os.path.exists(datapath):
            sys.exit('Could not find specified location for existing data: '+datapath)

    ###################################################################
    return knobs,dtstart,dtend,outputpath,extract,datapath,acdipole,OFB

############################################################################################################################################# END OF FUNCTION parse_options
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
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
#### OBTAIN KNOB DEFINITION USING THE ONLINE MODEL
def get_knob_circuits(knobs,dtstart,extract,outputpath,datapath):
    circuits_in_knob=[]
    OMoutputfile='knobs.madx'
    for i in range(len(knobs)):
        these_circuits=[]
        print ''
        print 'Obtaining circuits in knob:  '+knobs[i]
        print '------------------------------------------------'
        if extract==True:
            ##### CHECK FOR PRE-EXISTING OUTPUT, IF PRESENT MOVE IT TO AN 'OLD_' DIRECTORY NAME.
            if os.path.exists(outputpath+OMoutputfile):
                print 'Found prior output from Online Model'
                old_output_directory=outputpath+'old_output'
                if os.path.exists(old_output_directory):
                    pass
                else:
                    os.mkdir(old_output_directory)
                print 'Moving prior output to: '+old_output_directory
                print 'Even older extractions are being removed!!!'
                shutil.move(outputpath+OMoutputfile,old_output_directory+OMoutputfile)
            ##### RUN THE ONLINE MODEL FOR THE START OF THE SCAN
            print '---> Running online model'
            timestamp_OM = convert_to_OM_input_format(dtstart)
            #command='/afs/cern.ch/eng/lhc_online_model/dev/bin/lhc-model-extractor.sh -beam B1 -time \"'+timestamp_OM+'\" -Knames \"'+knobs[i]+'\" -Ktype KNOBVALUE -p '+outputpath
            command='/afs/cern.ch/eng/lhc_online_model/dev/bin/lhc-model-extractor.dev_knobs.sh -beam B1 -time \"'+timestamp_OM+'\" -Knames \"'+knobs[i]+'\" -Ktype KNOBVALUE -p '+outputpath
            #print command #useful for debugging
            os.system(command)
            #os.popen(command) # popen used to suppress output from online model command. makes things cleaner. can use os.system(command) if want to debug
        else:
            pass
        ########## READ THE NEWLY CREATED OUTPUT, OR IN-CASE NO EXTRACTION PERFORMED, CHECK THAT KNOB NAMES AND TIME-STAMPS MATCH UP TO THE EXPECTATION.
        if extract==True: # decide on path from extraction or previous output. Check the requested file actually exist.
            path=outputpath
        elif extract==False:
            path=datapath
        if not os.path.exists(path+OMoutputfile):
            print path+OMoutputfile
            sys.exit('Failed to find online model output')
            
        OMoutput=open(path+OMoutputfile,'r') ### READ THROUGH THE ONLINE MODEL OUTPUT FILE AND IDENTIFY BEAM-PROCESS, KNOB SETTING, AND THE CIRCUIT NAME;
        for line in OMoutput.readlines():
            if re.search('beamprocess',line):
                beam_process = line.partition('beamprocess')[2].partition('at ')[0].strip()
                print '---> Beam Process = '+beam_process
            elif re.search('! Extract Delta values for knob',line):
                thisknob = line.partition('values for knob')[2].partition('at time')[0].strip()
                if not thisknob==knobs[i]: # CHECK THE OUTPUT BEING EXAMINED IS ACTUALLY FOR THE REQUESTED KNOB
                    sys.exit('Knob name in pre-existing Online Model output does not match the requested knob.')
                filetime = line.partition('at time')[2].strip().strip('.').strip()
                print '---> Extraction at time: '+filetime
                dtfiletime = parse_timestamp(filetime)
                #if not dtstart==dtfiletime: # CHECK THE OUTPUT BEAING EXAMINED ACTUALLY CORRESPONDS TO THE REQUESTED START TIME
                #    sys.exit('Time stamps in pre-existing Online Model output does not match the requested start time.')
            elif re.search('! Value of knob',line):
                knob_setting = float(line.partition('=')[2].partition(';')[0].strip())
                knob_trim_var = line.partition('=')[0].strip() # get the variable name for the knob definition. will be used to check for circuits defined at zero in next step.
                print '---> Knob setting = '+str(knob_setting)
            elif re.search('Circuit name:',line):
                circuit_name = line.partition('Circuit name:')[2].strip().strip('.').strip('KICK').strip('/').strip()
                mad_var_name = line.partition('=')[0].strip()
                circuit_definition = line.partition(knob_trim_var)[0].partition('=')[2].partition(mad_var_name)[2].strip().strip('+').strip('-').strip('*').strip() ### Want to identify the knob definition for each circuit in order to exclude any which are defined at zero, e.g. in anglescan knobs.
                try:
                    circuit_definition=float(circuit_definition)
                except:
                    sys.exit('Problem identifying the circuit defininition in '+circuit_name)
                print '---> Found Circuit: '+circuit_name+'  with knob definition: '+str(circuit_definition)
                if circuit_definition==0.0:
                    print '---------------> Ommitting circuit: '+circuit_name+' from list of circuits as defined as zero in knob'
                    print '---------------> ('+line.strip()+')'
                else:
                    these_circuits.append(circuit_name)
        circuits_in_knob.append(these_circuits)
        OMoutput.close()

    return circuits_in_knob
###########################################################################################################################################################################
### TO UPDATE -> PANDAS
#### extract state of a list of circuits
def extract_circuit_state(list_of_circuits,dtstart,dtend): ## function to take the circuit_names found with online model, identify the corresponding timber variable and extract. Also make sure have an unambiguous determination of timber variable...
    print '---> Extracting data from timber...'
    ldb = pytimber.LoggingDB()
    temp_state_data=[]
    header=['#Time']
    for j in range(len(list_of_circuits)): ### this is the list of circuits corresponding to ONE of the input knobs.
        circuit_name=list_of_circuits[j]
        header.append(circuit_name)
        search_term='%'+circuit_name+'%STATE'  ### will search the timber database for varibles containing the circuit name and the STATE identifier
        search_result=ldb.search(search_term)
        if len(search_result)>1:
            sys.exit('STATE variable was not found unambiguously for search: '+search_term) ### making sure we only find a single timber variable which matches the search criteria. Don't know why this wouldn't be the case but probably as well to check.
        timber_variable=search_result[0] ### searching timber returns an array of matching variables. since have already determined there is only one match, take it explicitly out of the array
        print '---> Extracting '+timber_variable
        timber_st = convert_to_timber_input_format(dtstart)
        timber_et = convert_to_timber_input_format(dtend)
        this_logging_data = ldb.get(timber_variable,timber_st,timber_et,unixtime=False) ### HAVE RETRIEVED THE TIMBER DATA FOR THE REQUESTED CIRCUIT. UNIXTIME=FALSE RETURNS DATETIME OBJ IN LOCAL TIME (FOR NOW 03/2017, RdM says there is a UTC option, but unable to find it at present)

        ### this section is pretty unweildy and can definately be improved...

        for k in range(len(this_logging_data[timber_variable][0])):       ### Want to obtain a single list with all the trims of the circuits. complicated as while knob trims should all be synchronus, some circuits may be included in CO trims for example. 
            this_time=this_logging_data[timber_variable][0][k]            ### Want to merge into one big array, with common timestamps and settings for all circuits.
            this_state=this_logging_data[timber_variable][1][k]           ### Do this by appending all data for trims individucally, with blank entries for other circuits. then sorting the list on the datetime objects
            data_to_append=[this_time]                                    ### Then go through and merge common datetimes and pick out the non-common trims
            for l in range(len(list_of_circuits)):                   
                if l==j:
                    data_to_append.append(this_state)
                else:
                    data_to_append.append(None)
            temp_state_data.append(data_to_append)  
                                                     ### have data like    time, Circ1  , Circ2 , Circ3...
                                                     ###                   t1 ,  Idle , None , None  
                                                     ###                   t4 ,  Run  , None , None  
                                                     ###                   t1 ,  None , Idle , none
                                                     ###                   t2 ,  None , Run  , None  
                                                     ###                   t3 ,  None , Idle , None  
                                                     ###                   t4 ,  None , Run  , None
                                                     ###                   t1 ,  None , None , Idle
                                                     ###                   t4 ,  None , None , Run 


    sorted_state_data=sorted(temp_state_data,key=itemgetter(0)) ## sorting extracted data for all circuits by timestamp
                                                     ### have data like    time, Circ1  , Circ2 , Circ3...
                                                     ###                   t1 ,  Idle , None , None  -> first setting, all 3 circuits at idle 
                                                     ###                   t1 ,  None , Idle , none
                                                     ###                   t1 ,  None , None , Idle
                                                     ###                   t2 ,  None , Run  , None  -> None knob trim. eg Circ 2 is used for CO correction
                                                     ###                   t3 ,  None , Idle , None  -> Circ 2 CO correction finishes
                                                     ###                   t4 ,  Run  , None , None  -> Now have first real knob trim. All 3 circuits go to run, so have 3 time entries.
                                                     ###                   t4 ,  None , Run  , None
                                                     ###                   t4 ,  None , None , Run
#    sortedraw=sorted_state_data ### for debug
#    for z in range(len(sortedraw)):
#        print sortedraw[z]
#    sys.exit()

    for m in range(1,len(sorted_state_data)):
        if sorted_state_data[m][0]==sorted_state_data[m-1][0]:
            for n in range(len(sorted_state_data[m])):
                if sorted_state_data[m][n]==None:
                    sorted_state_data[m][n]=sorted_state_data[m-1][n]
                                                     ### have data like    time, Circ1  , Circ2 , Circ3...
                                                     ###                   t1 ,  Idle , None , None 
                                                     ###                   t1 ,  Idle , Idle , none
                                                     ###                   t1 ,  Idle , Idle , Idle
                                                     ###                   t2 ,  None , Run  , None 
                                                     ###                   t3 ,  None , Idle , None 
                                                     ###                   t4 ,  Run  , None , None 
                                                     ###                   t4 ,  Run  , Run  , None
                                                     ###                   t4 ,  Run  , Run  , Run
    state_data=[]
    state_data.append(header)
    for m in range(2,len(sorted_state_data)):
        if sorted_state_data[m][0]!=sorted_state_data[m-1][0]:
            state_data.append(sorted_state_data[m-1])
        if sorted_state_data[m][0]==sorted_state_data[m-1][0] and m==len(sorted_state_data):
            state_data.append(sorted_state_data[m])
                                                     ### have data like    time, Circ1  , Circ2 , Circ3...
                                                     ###                   t1 ,  Idle , Idle , Idle
                                                     ###                   t2 ,  None , Run  , None 
                                                     ###                   t3 ,  None , Idle , None 
                                                     ###                   t4 ,  Run  , Run  , Run    
    for m in range(2,len(state_data)):
        for n in range(len(state_data[m])):
            if state_data[m][n]==None:
                state_data[m][n]=state_data[m-1][n]
                                                     ### have data like    time, Circ1  , Circ2 , Circ3...
                                                     ###                   t1 ,  Idle , Idle , Idle
                                                     ###                   t2 ,  Idle , Run  , Idle 
                                                     ###                   t3 ,  Idle , Idle , Idle 
                                                     ###                   t4 ,  Run  , Run  , Run    


    return state_data

###########################################################################################################################################################################
def write_circuit_state(outputpath,knobname,this_state_data):
    knobname_for_file = convert_name(knobname)
    statefile=outputpath+'data.circuit_states.'+knobname_for_file+'.csv'
    print '---> Writing extracted state data to file:  '+statefile
    wout=open(statefile,'w')
    csvwout=csv.writer(wout,delimiter=',')
    for i in range(len(this_state_data)):
        if i==0:
            csvwout.writerow(this_state_data[i]) ### writing the header line to the output file
        else:
            #timestamp=this_state_data[i][0].strftime('%Y-%m-%d %H:%M:%S.%f')+'  ' ### output standard timestamp format rather than native datetime. neat formatting of this file will make checking for errors much easier
            timestamp=convert_to_data_output_format(this_state_data[i][0])+'  ' 
            datatowrite=[]
            datatowrite.append(timestamp)
            for j in range(1,len(this_state_data[i])):  ### adding all the state data to each row for writing out, but making up each csv entry to a standard length. improves readability of the output file for bug checking.
                somedata=this_state_data[i][j]
                for k in range(len(somedata),8):       
                    somedata=somedata+' '
                datatowrite.append(somedata)
            csvwout.writerow(datatowrite)
    wout.close()
    return
###########################################################################################################################################################################
def read_circuit_state(datapath,knobname):
    knobname_for_file = convert_name(knobname)
    statefile='data.circuit_states.'+knobname_for_file+'.csv'
    if os.path.exists(datapath+statefile):
        print '---> Retrieving existing state data from: '+datapath+statefile
    else:
        sys.exit('Could not find existing state data: '+datapath+statefile)
    ####### read the data extracted by previous run of the code, converting timestamps to datetime objects and removing whitespace included for formatting of the written file
    rdat=open(datapath+statefile,'r')
    csvrdat=csv.reader(rdat,delimiter=',',skipinitialspace=True)
    state_data=[]
    for row in csvrdat:
        if re.search('#',row[0]):
            data_to_append=[]
            for i in range(len(row)):
                data_to_append.append(row[i].strip())
        else:
            data_to_append=[]
            for i in range(len(row)):
                if i==0:
                    timestamp=row[i].strip()
                    #thisdt=dt.datetime.strptime(timestamp,'%Y-%m-%d %H:%M:%S.%f')
                    thisdt=convert_from_data_output_format(timestamp)
                    data_to_append.append(thisdt)
                else:
                    data_to_append.append(row[i].strip())
        state_data.append(data_to_append)
    rdat.close()
    if len(state_data)==0:
        sys.exit('Found state data file, but appears to be empty...')

    return state_data
###########################################################################################################################################################################
###########################################################################################################################################################################
def append_platteau(platteaustart,platteauend,platteaudat):
    thisplatteau=[]
    if platteaustart==None or platteauend==None:
        sys.exit('Something has gone very wrong with Platteau identification...')
    thisplatteau=[platteaustart,platteauend]
    platteaudat.append(thisplatteau)
    return platteaudat
############################################################################
def get_knob_platteaus(knobs,all_state_data,dtstart,dtend):  ### Aim here is to return timestamps of the knob trims.
    knob_platteaus=[]                                        ### Without a better way, will identify a knob trim as a common change of state of all circuits in the knob
    for i in range(len(knobs)):                              ### This should exclude most cases where knob circuits are included in an orbit feed-back
        knobname=knobs[i]                                    ### Should also include a check of the knob setting, to avoid a case where all knob circuits are used in the OFB...
        print '\n\n'                                         ### what to do if there is a constant OFB with all knobs included...? - should be caught by the previous check... 
        print 'Checking trim platteaus for knob: '+knobname
        print '----------------------------------------------------------------------'
        this_knob_state_data = all_state_data[i]  ### selecting the state data corresponding to the knob being looked at.
        this_knob_common_state_changes=[]
        for j in range(1,len(this_knob_state_data)):
            states=set(this_knob_state_data[j][1:])  ### take a list of all the circuit states at a given time. Convert to a set. As set removes duplicates, if the length of the new set is 1 all circuits have a common state change at that time.
            if len(states)==1:
                this_knob_common_state_changes.append([this_knob_state_data[j][0],list(states)[0]]) ## for common trims append time and trim state

        this_knob_platteaus=[]
        platteaustart=dtstart
        currentstate='IDLE'
        for j in range(len(this_knob_common_state_changes)):
            newstate=this_knob_common_state_changes[j][1]
            if newstate!=currentstate:
                if currentstate=='RUNNING' and newstate=='IDLE':
                    platteaustart=this_knob_common_state_changes[j][0]
                    if j==len(this_knob_common_state_changes)-1:
                        platteauend=dtend
                        this_knob_platteaus = append_platteau(platteaustart,platteauend,this_knob_platteaus)
                elif currentstate=='IDLE' or currentstate=='ARMED':
                    if newstate=='RUNNING':
                        platteauend=this_knob_common_state_changes[j][0]
                        this_knob_platteaus = append_platteau(platteaustart,platteauend,this_knob_platteaus)
                        platteaustart=None
                        platteauend=None
                currentstate=newstate
        
        print '---> Found '+str(len(this_knob_platteaus))+' knob trim platteaus'
        for j in range(len(this_knob_platteaus)):
            starttime=convert_to_data_output_format(this_knob_platteaus[j][0])
            endtime=convert_to_data_output_format(this_knob_platteaus[j][1])
            print '---> '+starttime+' -> '+endtime

        knob_platteaus.append(this_knob_platteaus)
    return knob_platteaus

###########################################################################################################################################################################
def progressbar(it, prefix="", size=60):
    count = len(it)
    def _show(_i):
        x = int(size*_i/count)
        sys.stdout.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), _i, count))
        sys.stdout.flush()

    _show(0)
    for i, item in enumerate(it):
        yield item
        _show(i+1)
    sys.stdout.write("\n")
    sys.stdout.flush()
###########################################################################################################################################################################            
def get_mcb_states(knobs,dtstart,dtend):
    all_mcb_states=[]
    ldb = pytimber.LoggingDB()
    for i in range(len(knobs)):
        print 'KNOB = '+knobs[i]
        knob_type=knobs[i].partition('/')[0]
        if knob_type=='LHCBEAM':
            print '---> knob type identified as LHCBEAM. Will check for MCBX; MCB.B1; MCB.B2'
            search_terms=['%RCBX%STATE','%RCB%B1%STATE','%RCB%B2%STATE']          
        elif knob_type=='LHCBEAM1':
            print '---> knob type identified as LHCBEAM. Will check for MCBX; MCB.B1'
            search_terms=['%RCBX%STATE','%RCB%.B1%STATE']
        elif knob_type=='LHCBEAM2':
            print '---> knob type identified as LHCBEAM. Will check for MCBX; MCB.B2'
            search_terms=['%RCBX%STATE','%RCB%.B2%STATE']     

        circuits=[]
        for j in range(len(search_terms)):
            search_result=ldb.search(search_terms[j])
            for k in range(len(search_result)):
                circuits.append(search_result[k])
        print '---> Found '+str(len(circuits))+' relevant MCB circuits'
        
        if len(circuits)==0:
            sys.exit('FAILED TO FIND MCB??')
        
        timber_st = convert_to_timber_input_format(dtstart)
        timber_et = convert_to_timber_input_format(dtend)
        
        mcb_states=[]
        for j in progressbar(range(len(circuits)), "---> Extracting circuit states: ", 100):
            circuit_name=circuits[j]
            this_logging_data = ldb.get(circuit_name,timber_st,timber_et,unixtime=False)
            #var=ldb.search("%RCBYHS4.R5B1%STATE%")
            #var=var[0]
            #this_logging_data = ldb.get(var,timber_st,timber_et,unixtime=False)
            #time=this_logging_data[var][0]
            timedat=this_logging_data[circuit_name][0]
            statedat=this_logging_data[circuit_name][1]
            for k in range(len(timedat)):
                mcbdat=[timedat[k],statedat[k],circuit_name]
                mcb_states.append(mcbdat)
#            if j==20:
#                break
#####        mcb_states=sorted(mcb_states,key=itemgetter(0))
    all_mcb_states.append(mcb_states)
    return all_mcb_states
###########################################################################################################################################################################
###########################################################################################################################################################################
def write_mcb_states(outputpath,knobs,all_mcb_states):
    for i in range(len(knobs)):
        knobname_for_file = convert_name(knobs[i])
        mcbfile=outputpath+'data.all_mcb_states.'+knobname_for_file+'.csv'
        print '---> Writing mcb states to file:  '+mcbfile
        wout=open(mcbfile,'w')
        csvwout=csv.writer(wout,delimiter=',')
        this_mcb_states=all_mcb_states[i] ## selecting the states array corresponding to the knob in question
        for j in range(len(this_mcb_states)):
            timestamp=convert_to_data_output_format(this_mcb_states[j][0])
            state=this_mcb_states[j][1]
            circuit=this_mcb_states[j][2]
            datatowrite=[timestamp,state,circuit]
            csvwout.writerow(datatowrite)
        wout.close()
        
        mcbrunfile=outputpath+'data.mcb.runfile.'+knobname_for_file+'.csv'
        wout=open(mcbrunfile,'w')
        csvwout=csv.writer(wout,delimiter=',')
        tempmcbstate=cp.deepcopy(this_mcb_states)
        sorted_mcb_states=sorted(tempmcbstate,key=itemgetter(0))
        for j in range(len(sorted_mcb_states)):
            state=this_mcb_states[j][1]
            if state=='RUNNING':
                timestamp=convert_to_data_output_format(this_mcb_states[j][0])
                circuit=this_mcb_states[j][2]
                datatowrite=[timestamp,state,circuit]
                csvwout.writerow(datatowrite)
        wout.close()

    return
###########################################################################################################################################################################
def read_mcb_states(datapath,knobs):
    all_mcb_states=[]
    for i in range(len(knobs)):
        print 'KNOB = '+knobs[i]
        mcb_states=[]
        knobname_for_file = convert_name(knobs[i])
        print knobname_for_file
        print datapath+'data.all_mcb_states.'
        mcbfile=datapath+'data.all_mcb_states.'+knobname_for_file+'.csv'
        if os.path.exists(mcbfile):
            print '---> Reading mcb states from file:  '+mcbfile
        else:
            sys.exit('Could not find existing mcb state data: '+mcbfile)
        rdat=open(mcbfile)
        csvrdat=csv.reader(rdat,delimiter=',',skipinitialspace=True)
        for row in csvrdat:
            if re.search('#',row[0]):
                continue
            timestamp=row[0].strip()
            thisdt=convert_from_data_output_format(timestamp)
            state=row[1].strip()
            circuit=row[2].strip()
            mcb_states.append([thisdt,state,circuit])
        rdat.close()
        if len(mcb_states)==0:
            sys.exit('Found data file, but appears to be empty...')
        all_mcb_states.append(mcb_states)
    return all_mcb_states
###########################################################################################################################################################################
def is_knob_trim(knob_platteaus,mcb_trim_time):
    for k in range(len(knob_platteaus)):
        knob_plat_end=knob_platteaus[k][1]
        if mcb_trim_time==knob_plat_end:
            return True
    return False
###########################################################################################
def is_already_platteau_start(data_platteaus,mcb_trim_end_time):
    for k in range(len(data_platteaus)):
        data_plat_start=data_platteaus[k][0]
        if mcb_trim_end_time==data_plat_start:
            return True
    return False
###########################################################################################
def update_data_platteaus(data_platteaus,mcb_trim_end_time):
    for k in range(len(data_platteaus)):
        data_plat_start=data_platteaus[k][0]
        data_plat_end=data_platteaus[k][1]
        if mcb_trim_end_time>data_plat_start and mcb_trim_end_time<data_plat_end:
            print '---> Identified CO correction at '+convert_to_data_output_format(mcb_trim_end_time)
            data_platteaus[k][0]=mcb_trim_end_time ### update start time of the data platteau to account for the CO trim
    return data_platteaus
###########################################################################################
def check_for_CO_correction(knobs,all_knob_platteaus,all_mcb_states,dtstart,dtend):
    all_data_platteaus=[]
    for i in range(len(knobs)):
        print 'KNOB = '+knobs[i]
        mcb_states=all_mcb_states[i]
        knob_platteaus=all_knob_platteaus[i]
        data_platteaus=cp.deepcopy(knob_platteaus) ### Mirror the knob platteaus. this is the list we will edit when find CO corrections.
        for j in range(len(mcb_states)):
            mcb_trim_time=mcb_states[j][0]
            state=mcb_states[j][1]
            circuit=mcb_states[j][2]
            if state!='RUNNING': ### only looking for mcb states which are running as these indicate the start of a trim
                continue
            if is_knob_trim(knob_platteaus,mcb_trim_time)==True: ### if circuit trim corresponds to the knob trim, want to ignore it
                continue
            next_mcb_time=mcb_states[j+1][0]
            next_mcb_state=mcb_states[j+1][1]
            next_mcb_circuit=mcb_states[j+1][2]
            if next_mcb_circuit!=circuit:
                sys.exit('Final trim of circuit: '+circuit+' within the specified measurement period is to \' '+state+' \'. Occurs at '+mcb_trim_time+'. Probably means have a bad end time selected for the scan.')
            if next_mcb_state!='IDLE':
                sys.exit('MCB state following \'Running\' is not \'IDLE\'. Why??')
            mcb_trim_end_time=mcb_states[j+1][0]
            if is_already_platteau_start(data_platteaus,mcb_trim_end_time)==True: ### if circuit trim corresponds to a time which is already the start of a data platteau (data platteau, initially this is the knob trim start, but will already be updated for earlier found trims), want to ignore it
                continue
            else:
                data_platteaus = update_data_platteaus(data_platteaus,mcb_trim_end_time)
        if len(data_platteaus)!=len(knob_platteaus):
            sys.exit('CO correction search has changed number of Platteaus.... ')

        print 'New data platteaus for knob = '+knobs[i]
        for j in range(len(data_platteaus)):
            originalstarttime=convert_to_data_output_format(knob_platteaus[j][0])
            starttime=convert_to_data_output_format(data_platteaus[j][0])
            endtime=convert_to_data_output_format(data_platteaus[j][1])
            if data_platteaus[j][1]!=knob_platteaus[j][1]:
                sys.exit('Seem to have changed end time of a platteau. Shouldn\'t happen.')
            print '---> '+originalstarttime+'  --->  '+starttime+' -> '+endtime

        all_data_platteaus.append(data_platteaus)
    
    return all_data_platteaus
###########################################################################################################################################################################
def write_platteaus(outputpath,knobs,all_knob_platteaus,all_data_platteaus):
    for i in range(len(knobs)):
        knobname_for_file = convert_name(knobs[i])
        platfile=outputpath+'data.platteaus.'+knobname_for_file+'.csv'
        print '---> Writing platteaus to file:  '+platfile
        wout=open(platfile,'w')
        csvwout=csv.writer(wout,delimiter=',')
        header=['#knob_plat','data_plat_start','data_plat_end']
        csvwout.writerow(header)
        this_knob_plat=all_knob_platteaus[i] ## selecting the states array corresponding to the knob in question
        this_data_plat=all_data_platteaus[i]
        for j in range(len(this_data_plat)):
            knob_plat_start=convert_to_data_output_format(this_knob_plat[j][0])
            data_plat_start=convert_to_data_output_format(this_data_plat[j][0])
            data_plat_end  =convert_to_data_output_format(this_data_plat[j][1])
            datatowrite=[knob_plat_start,data_plat_start,data_plat_end]
            csvwout.writerow(datatowrite)
    wout.close()
    return
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################


###########################################################################################################################################################################
#### MAIN FUNCTION FOR IRNL SCAN ##################################################################################################### START OF FUNCTION main_get_platteaus
###########################################################################################################################################################################
def main_get_platteaus(knobs,dtstart,dtend,extract,outputpath,datapath,OFB):
    ###################################################################
    circuits_in_knob = get_knob_circuits(knobs,dtstart,extract,outputpath,datapath) ## use online model to check which circuits are included in the knobs
    ###################################################################
    all_state_data=[]
    for i in range(len(circuits_in_knob)):
        print '\n\n'
        print 'Obtaining circuit state data for knob '+knobs[i]
        print '-------------------------------------------------------------'
        if extract==True:
            state_data = extract_circuit_state(circuits_in_knob[i],dtstart,dtend)
            write_circuit_state(outputpath,knobs[i],state_data)
        elif extract==False:
            state_data = read_circuit_state(datapath,knobs[i])
        all_state_data.append(state_data)
    ###################################################################
    all_knob_platteaus = get_knob_platteaus(knobs,all_state_data,dtstart,dtend)
    ###################################################################
    if OFB==True:  ### If automatic orbit feed-back is operational don't want to check for individual trims of corrector circuits, as they will alway be acting to keep orbit good.
        all_data_platteaus=cp.deepcopy(all_knob_platteaus)
    else:    ### However if there is no automatic OFB, need to check for any potential manual orbit correction which could have been performed during one of the knob-trim platteaus
        print '\n\n'
        print 'Obtaining state data for all orbit corrector circuits in the machine'
        print '--------------------------------------------------------------------' #################### DEBUG
        
        if extract==True:
            all_mcb_states = get_mcb_states(knobs,dtstart,dtend)
            write_mcb_states(outputpath,knobs,all_mcb_states)
        else:
            all_mcb_states = read_mcb_states(datapath,knobs)
        '''
        if extract==True:
            all_mcb_states = read_mcb_states(outputpath,knobs)                         #################### DEBUG
        else:
            all_mcb_states = read_mcb_states(datapath,knobs)                           #################### DEBUG  
        '''
        print '\n\n'
        print 'Checking for individual trims of MCB within knob platteaus (e.g. manaul orbit corrections)'
        print '-------------------------------------------------------------------------------------------'
        all_data_platteaus = check_for_CO_correction(knobs,all_knob_platteaus,all_mcb_states,dtstart,dtend) ### New set of platteaus, where any manual CO corr results in earlier data in the platteau being excluded
    write_platteaus(outputpath,knobs,all_knob_platteaus,all_data_platteaus)
    ###################################################################
    return circuits_in_knob,all_data_platteaus

###########################################################################################################################################################################
#### END OF FUNCTION main ############################################################################################################## END OF FUNCTION main_get_platteaus
###########################################################################################################################################################################





###########################################################################################################################################################################
#### MAIN FUNCTION FOR IRNL SCAN ###################################################################################################### START OF FUNCTION main_analysis_bbq
###########################################################################################################################################################################
def main_analysis_bbq():
    ###################################################################
    return
###########################################################################################################################################################################
#### END OF FUNCTION main ############################################################################################################### END OF FUNCTION main_analysis_bbq
###########################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
###########################################################################################################################################################################
#### MAIN FUNCTION FOR IRNL SCAN ###################################################################################################### START OF FUNCTION main_analysis_acd
###########################################################################################################################################################################
def main_analysis_acd():
    ###################################################################
    return
###########################################################################################################################################################################
#### END OF FUNCTION main ############################################################################################################### END OF FUNCTION main_analysis_acd
###########################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
###########################################################################################################################################################################
#### MAIN FUNCTION FOR IRNL SCAN ############################################################################################################# START OF FUNCTION main_model
###########################################################################################################################################################################
def main_model():
    ###################################################################
    return
###########################################################################################################################################################################
#### END OF FUNCTION main ###################################################################################################################### END OF FUNCTION main_model
###########################################################################################################################################################################


######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
def get_orbit_meta(orbitvar,dtstart,dtend,ldb):
    print '\nExtracting orbit metadata for '+orbitvar
    orbit_meta=ldb.getMetaData(orbitvar)
    orbit_meta_times=orbit_meta[orbitvar][0]
    orbit_meta_bpms=orbit_meta[orbitvar][1]
    for j in range(len(orbit_meta_times)):
        time_of_BPM_list=unix_to_dtime_converter(orbit_meta_times[j])
        print '---> Found metadata from: '+convert_to_data_output_format(time_of_BPM_list)
        if j==len(orbit_meta_times)-1:
            print 'Using most recent BPM metadata, from '+convert_to_data_output_format(time_of_BPM_list)
            meta_to_use=orbit_meta_bpms[j]
            continue
        next_time_of_BPM_list=unix_to_dtime_converter(orbit_meta_times[j+1])
        if time_of_BPM_list<dtstart and next_time_of_BPM_list>dtstart:
            if dtstart<next_time_of_BPM_list and dtend>next_time_of_BPM_list:
                sys.exit('BPM metatdata seems to change during the scan??? Think something probably went wrong here...')
            else:
                print '---> Using BPM metadata from '+convert_to_data_output_format(time_of_BPM_list)
                print '     (next metadata is from '+convert_to_data_output_format(next_time_of_BPM_list)+')'
                meta_to_use=orbit_meta_bpms[j]
                for k in range(j,len(orbit_meta_times)):
                    print '---> Found metadata from: '+convert_to_data_output_format(orbit_meta_times[k])
                continue
    return meta_to_use
##########################################################################################################
def get_orbit_data(orbitvar,dtstart,dtend,ldb):
    print '\nExtracting orbit data for '+orbitvar
    
    timber_st = convert_to_timber_input_format(dtstart)
    #dtime=dt.timedelta(0,10,0)
    #timber_et = convert_to_timber_input_format(dtstart+dtime)
    timber_et = convert_to_timber_input_format(dtend)

    timber_orbit_data=ldb.get(orbitvar,timber_st,timber_et,unixtime=False)
    timber_orbit_data_times=timber_orbit_data[orbitvar][0]
    timber_orbit_data_pos=timber_orbit_data[orbitvar][1]
    orbit_data=[]
    for j in range(len(timber_orbit_data_times)):
        datarow=[]
        time=timber_orbit_data_times[j]
        datarow.append(time)
        for k in range(len(timber_orbit_data_pos[j])):
            datarow.append(timber_orbit_data_pos[j][k])
        orbit_data.append(datarow)
    return orbit_data
##########################################################################################################
def write_orbit_data(outputpath,orbitvar,meta_to_use,orbit_data):
    varforfile=convert_name(orbitvar)
    filename=outputpath+'data.orbit.raw.'+varforfile+'.csv'
    print '---> Writing raw orbit data to file:  '+filename
    wout=open(filename,'w')
    csvwout=csv.writer(wout,delimiter=',')
    header=[]
    header.append('#timestamp(data in um)')
    for j in range(len(meta_to_use)):
        header.append(meta_to_use[j])
    csvwout.writerow(header)
    for j in range(len(orbit_data)):
        if len(header)!=len(orbit_data[j]):
            print len(header)
            print len(orbit_data[j])
            sys.exit('Undefined BPMs?')
        else:
            datarow=[]
            timestamp=convert_to_data_output_format(orbit_data[j][0])
            datarow.append(timestamp)
            for k in range(1,len(orbit_data[j])):
                datarow.append(float(orbit_data[j][k]))
            if len(header)!=len(datarow):
                sys.exit('Undefined BPMs.2?')
            csvwout.writerow(datarow)
    wout.close()
    return
##########################################################################################################
def extract_orbit(outputpath,orbitvars,dtstart,dtend,ldb):
    all_meta_to_use=[]
    all_orbit_data=[]
    for i in range(len(orbitvars)):
        orbitvar=orbitvars[i]
        meta_to_use = get_orbit_meta(orbitvar,dtstart,dtend,ldb)
        orbit_data = get_orbit_data(orbitvar,dtstart,dtend,ldb)
        write_orbit_data(outputpath,orbitvar,meta_to_use,orbit_data)
        all_meta_to_use.append(meta_to_use)
        all_orbit_data.append(orbit_data)
    return all_meta_to_use,all_orbit_data
##########################################################################################################
def define_orbit_vars(): ### IF REDEFINE THIS MUST ALSO REDEFINE THE FUNCTION BELOW 'merge_orbit_properties()'
    orbitHvar='LHC.BOFSU:POSITIONS_H'
    orbitVvar='LHC.BOFSU:POSITIONS_V'
    orbitvars=[orbitHvar,orbitVvar]
    return orbitvars
##########################################################################################################
def merge_orbit_properties(all_orbit_properties):
    Hdata=cp.deepcopy(all_orbit_properties[0])
    Hheader=Hdata[0]
    del Hdata[0]
    Vdata=cp.deepcopy(all_orbit_properties[1])
    Vheader=Vdata[0]
    del Vdata[0]

    main_data=[]
    detailed_data=[]

    main_header=['#time','Hmean_arc_orbit_b1','Vmean_arc_orbit_b1','Hrms_arc_orbit_b1','Vrms_arc_orbit_b1','Hxing_ir1b1','Vxing_ir1b1','Hxing_ir5b1','Vxing_ir5b1','Hmean_arc_orbit_b2','Vmean_arc_orbit_b2','Hrms_arc_orbit_b2','Vrms_arc_orbit_b2','Hxing_ir1b2','Vxing_ir1b2','Hxing_ir5b2','Vxing_ir5b2']
    detailed_header=['#time','Horb1l1b1','Horb1r1b1','Horb2l1b1','Horb2r1b1','Vorb1l1b1','Vorb1r1b1','Vorb2l1b1','Vorb2r1b1','Horb1l5b1','Horb1r5b1','Horb2l5b1','Horb2r5b1','Vorb1l5b1','Vorb1r5b1','Vorb2l5b1','Vorb2r5b1','Horb1l1b2','Horb1r1b2','Horb2l1b2','Horb2r1b2','Vorb1l1b2','Vorb1r1b2','Vorb2l1b2','Vorb2r1b2','Horb1l5b2','Horb1r5b2','Horb2l5b2','Horb2r5b2','Vorb1l5b2','Vorb1r5b2','Vorb2l5b2','Vorb2r5b2']

    for h in range(len(Hdata)):
        Htime=Hdata[h][0]
        Hmean_arc_orbit_b1=Hdata[h][1]
        Hmean_arc_orbit_b2=Hdata[h][2]
        Hrms_arc_orbit_b1=Hdata[h][3]
        Hrms_arc_orbit_b2=Hdata[h][4]
        Hxing_ir1b1=Hdata[h][5]
        Hxing_ir1b2=Hdata[h][6]
        Hxing_ir5b1=Hdata[h][7]
        Hxing_ir5b2=Hdata[h][8]
        Horb1l1b1=Hdata[h][9]
        Horb1r1b1=Hdata[h][10]
        Horb2l1b1=Hdata[h][11]
        Horb2r1b1=Hdata[h][12]
        Horb1l1b2=Hdata[h][13]
        Horb1r1b2=Hdata[h][14]
        Horb2l1b2=Hdata[h][15]
        Horb2r1b2=Hdata[h][16]
        Horb1l5b1=Hdata[h][17]
        Horb1r5b1=Hdata[h][18]
        Horb2l5b1=Hdata[h][19]
        Horb2r5b1=Hdata[h][20]
        Horb1l5b2=Hdata[h][21]
        Horb1r5b2=Hdata[h][22]
        Horb2l5b2=Hdata[h][23]
        Horb2r5b2=Hdata[h][24]
        found_matching_Vdata=False
        for v in range(len(Vdata)):
            Vtime=Vdata[v][0]
            if Htime==Vtime:
                found_matching_Vdata=True
                Vmean_arc_orbit_b1=Vdata[v][1]
                Vmean_arc_orbit_b2=Vdata[v][2]
                Vrms_arc_orbit_b1=Vdata[v][3]
                Vrms_arc_orbit_b2=Vdata[v][4]
                Vxing_ir1b1=Vdata[v][5]
                Vxing_ir1b2=Vdata[v][6]
                Vxing_ir5b1=Vdata[v][7]
                Vxing_ir5b2=Vdata[v][8]
                Vorb1l1b1=Vdata[v][9]
                Vorb1r1b1=Vdata[v][10]
                Vorb2l1b1=Vdata[v][11]
                Vorb2r1b1=Vdata[v][12]
                Vorb1l1b2=Vdata[v][13]
                Vorb1r1b2=Vdata[v][14]
                Vorb2l1b2=Vdata[v][15]
                Vorb2r1b2=Vdata[v][16]
                Vorb1l5b1=Vdata[v][17]
                Vorb1r5b1=Vdata[v][18]
                Vorb2l5b1=Vdata[v][19]
                Vorb2r5b1=Vdata[v][20]
                Vorb1l5b2=Vdata[v][21]
                Vorb1r5b2=Vdata[v][22]
                Vorb2l5b2=Vdata[v][23]
                Vorb2r5b2=Vdata[v][24]
                del Vdata[v]
                break
        if found_matching_Vdata==False:
            print 'nooooooooooooooooooo'
            Vmean_arc_orbit_b1=float('NaN')
            Vmean_arc_orbit_b2=float('NaN')
            Vrms_arc_orbit_b1=float('NaN')
            Vrms_arc_orbit_b2=float('NaN')
            Vxing_ir1b1=float('NaN')
            Vxing_ir1b2=float('NaN')
            Vxing_ir5b1=float('NaN')
            Vxing_ir5b2=float('NaN')
            Vorb1l1b1=float('NaN')
            Vorb1r1b1=float('NaN')
            Vorb2l1b1=float('NaN')
            Vorb2r1b1=float('NaN')
            Vorb1l1b2=float('NaN')
            Vorb1r1b2=float('NaN')
            Vorb2l1b2=float('NaN')
            Vorb2r1b2=float('NaN')
            Vorb1l5b1=float('NaN')
            Vorb1r5b1=float('NaN')
            Vorb2l5b1=float('NaN')
            Vorb2r5b1=float('NaN')
            Vorb1l5b2=float('NaN')
            Vorb1r5b2=float('NaN')
            Vorb2l5b2=float('NaN')
            Vorb2r5b2=float('NaN')

        thismaindata=[Htime,Hmean_arc_orbit_b1,Vmean_arc_orbit_b1,Hrms_arc_orbit_b1,Vrms_arc_orbit_b1,Hxing_ir1b1,Vxing_ir1b1,Hxing_ir5b1,Vxing_ir5b1,Hmean_arc_orbit_b2,Vmean_arc_orbit_b2,Hrms_arc_orbit_b2,Vrms_arc_orbit_b2,Hxing_ir1b2,Vxing_ir1b2,Hxing_ir5b2,Vxing_ir5b2]
        thisdetaileddata=[Htime,Horb1l1b1,Horb1r1b1,Horb2l1b1,Horb2r1b1,Vorb1l1b1,Vorb1r1b1,Vorb2l1b1,Vorb2r1b1,Horb1l5b1,Horb1r5b1,Horb2l5b1,Horb2r5b1,Vorb1l5b1,Vorb1r5b1,Vorb2l5b1,Vorb2r5b1,Horb1l1b2,Horb1r1b2,Horb2l1b2,Horb2r1b2,Vorb1l1b2,Vorb1r1b2,Vorb2l1b2,Vorb2r1b2,Horb1l5b2,Horb1r5b2,Horb2l5b2,Horb2r5b2,Vorb1l5b2,Vorb1r5b2,Vorb2l5b2,Vorb2r5b2]
        main_data.append(thismaindata)
        detailed_data.append(thisdetaileddata)

    for v in range(len(Vdata)): ### have already gone through and deleted any entries from Vdata that match up with H data. 
                                ### Hdata entries with no parallel in Vdata set all V parameters to nan and append. 
                                ### Now want to go through an remaining Vdata entries which don't have a parallel in Hdata and do the same. 
                                ### Only now know they have alredy matched up all common entries so don't need to worry about looping through the Hdata.
        Vtime=Vdata[v][0]
        Vmean_arc_orbit_b1=Vdata[v][1]
        Vmean_arc_orbit_b2=Vdata[v][2]
        Vrms_arc_orbit_b1=Vdata[v][3]
        Vrms_arc_orbit_b2=Vdata[v][4]
        Vxing_ir1b1=Vdata[v][5]
        Vxing_ir1b2=Vdata[v][6]
        Vxing_ir5b1=Vdata[v][7]
        Vxing_ir5b2=Vdata[v][8]
        Vorb1l1b1=Vdata[v][9]
        Vorb1r1b1=Vdata[v][10]
        Vorb2l1b1=Vdata[v][11]
        Vorb2r1b1=Vdata[v][12]
        Vorb1l1b2=Vdata[v][13]
        Vorb1r1b2=Vdata[v][14]
        Vorb2l1b2=Vdata[v][15]
        Vorb2r1b2=Vdata[v][16]
        Vorb1l5b1=Vdata[v][17]
        Vorb1r5b1=Vdata[v][18]
        Vorb2l5b1=Vdata[v][19]
        Vorb2r5b1=Vdata[v][20]
        Vorb1l5b2=Vdata[v][21]
        Vorb1r5b2=Vdata[v][22]
        Vorb2l5b2=Vdata[v][23]
        Vorb2r5b2=Vdata[v][24]
        
        Hmean_arc_orbit_b1=float('NaN')
        Hmean_arc_orbit_b2=float('NaN')
        Hrms_arc_orbit_b1=float('NaN')
        Hrms_arc_orbit_b2=float('NaN')
        Hxing_ir1b1=float('NaN')
        Hxing_ir1b2=float('NaN')
        Hxing_ir5b1=float('NaN')
        Hxing_ir5b2=float('NaN')
        Horb1l1b1=float('NaN')
        Horb1r1b1=float('NaN')
        Horb2l1b1=float('NaN')
        Horb2r1b1=float('NaN')
        Horb1l1b2=float('NaN')
        Horb1r1b2=float('NaN')
        Horb2l1b2=float('NaN')
        Horb2r1b2=float('NaN')
        Horb1l5b1=float('NaN')
        Horb1r5b1=float('NaN')
        Horb2l5b1=float('NaN')
        Horb2r5b1=float('NaN')
        Horb1l5b2=float('NaN')
        Horb1r5b2=float('NaN')
        Horb2l5b2=float('NaN')
        Horb2r5b2=float('NaN')

        thismaindata=[Htime,Hmean_arc_orbit_b1,Vmean_arc_orbit_b1,Hrms_arc_orbit_b1,Vrms_arc_orbit_b1,Hxing_ir1b1,Vxing_ir1b1,Hxing_ir5b1,Vxing_ir5b1,Hmean_arc_orbit_b2,Vmean_arc_orbit_b2,Hrms_arc_orbit_b2,Vrms_arc_orbit_b2,Hxing_ir1b2,Vxing_ir1b2,Hxing_ir5b2,Vxing_ir5b2]
        thisdetaileddata=[Htime,Horb1l1b1,Horb1r1b1,Horb2l1b1,Horb2r1b1,Vorb1l1b1,Vorb1r1b1,Vorb2l1b1,Vorb2r1b1,Horb1l5b1,Horb1r5b1,Horb2l5b1,Horb2r5b1,Vorb1l5b1,Vorb1r5b1,Vorb2l5b1,Vorb2r5b1,Horb1l1b2,Horb1r1b2,Horb2l1b2,Horb2r1b2,Vorb1l1b2,Vorb1r1b2,Vorb2l1b2,Vorb2r1b2,Horb1l5b2,Horb1r5b2,Horb2l5b2,Horb2r5b2,Vorb1l5b2,Vorb1r5b2,Vorb2l5b2,Vorb2r5b2]
        main_data.append(thismaindata)
        detailed_data.append(thisdetaileddata)

    ### Now have list of all H and V orbit data. Where one doesn't have a companion in the other set of data have filled with 'NaN'. But H and V data is unsorted. Want to sort the array now by timestamp.

    main_data=sorted(main_data,key=itemgetter(0))
    detailed_data=sorted(detailed_data,key=itemgetter(0))

    main_data.insert(0,main_header)
    detailed_data.insert(0,detailed_header)

    return main_data,detailed_data

##########################################################################################################
def read_orbit(datapath,orbitvars):
    all_meta_to_use=[]
    all_orbit_data=[]
    for i in range(len(orbitvars)):
        orbitvar=orbitvars[i]
        orbit_data=[]
        varforfile=convert_name(orbitvar)
        filename=datapath+'data.orbit.raw.'+varforfile+'.csv'
        if os.path.exists(filename):
            print '---> Reading raw orbit data from file:  '+filename
        else:
            sys.exit('Could not find existing orbit data: '+filename)
        rdat=open(filename,'r')
        csvrdat=csv.reader(rdat,delimiter=',',skipinitialspace=True)
        for row in csvrdat:
            if re.search('#timestamp',row[0]):
                meta_to_use=row[1:]
                continue
            else:
                datarow=[]
                timestamp=row[0]
                dtobj=convert_from_data_output_format(timestamp)
                datarow.append(dtobj)
                for z in range(1,len(row)):
                    datarow.append(float(row[z]))
            orbit_data.append(datarow)
        rdat.close()
        all_meta_to_use.append(meta_to_use)
        all_orbit_data.append(orbit_data)
    return all_meta_to_use,all_orbit_data
##########################################################################################################
##########################################################################################################
def is_arc_BPM(name):
    label,ref,beam=name.split('.')
    if label!='BPM':
        return 0
    ref=re.sub('[BA]','',ref)   # Some BPMs have a B or A in the ref: BPMWA.B5L4.B1, BPMWA.A5L4.B1
    ref=int(re.split('[RL]',ref)[0])
    try:
        if ref>10:
            return True
        else:
            return False
    except:
        print '-------> Failed to identify if '+name+' is an ARC BPM. Will proceed on assumption it is not.'
        return False
##########################################################################################################
def get_inferred_orbit(orbitvars,all_meta_to_use,all_orbit_data):
    all_orbit_properties=[]
    for p in range(len(orbitvars)): ### p defines the plane in accordance with the orbit variables defined in function define_orbit_vars()
       print '---> Calcultating orbit properties vs time from data '+orbitvars[p]
       meta_to_use=all_meta_to_use[p]
       orbit_data=all_orbit_data[p]
       orbit_properties=[]
       header=['#time','mean_arc_orbit_b1','mean_arc_oribt_b2','rms_arc_orbit_b1','rms_arc_orbit_b2','xing_ir1b1','xing_ir1b2','xing_ir5b1','xing_ir5b2','orb1l1b1','orb1r1b1','orb2l1b1','orb2r1b1','orb1l1b2','orb1r1b2','orb2l1b2','orb2r1b2','orb1l5b1','orb1r5b1','orb2l5b1','orb2r5b1','orb1l5b2','orb1r5b2','orb2l5b2','orb2r5b2']
       orbit_properties.append(header)

       b1_arc_nan=False
       b2_arc_nan=False
              
       for t in progressbar(range(len(orbit_data)), "------> Time stamps remaining: ", 100):
           time=orbit_data[t][0]
           #### lalala
           arc_orbit_b1=[]
           arc_orbit_b2=[]
           for b in range(len(meta_to_use)): ### loop through all the bpms
               this_bpm=meta_to_use[b]
               only_orbit_data=orbit_data[t][1:] ### strip out the timestamp entry in the row so that entry 'b' in meta_to_use, corresponds to orbit 'b' in the orbit_data list. otherwise have offset of 1 due to extra timestamp entry in the orbit data rows.
               try:
                   this_orbit=float(only_orbit_data[b])/(1.0e3) ## read bpm data, and convert from um to mm. set to nan if fail to get out an orbit value
               except:
                   this_orbit=float('NaN')
               ################################################ ### obtaining the orbit at some of the interesting bpms in the irs
               if this_bpm=='BPMSW.1L5.B1':
                   orb1l5b1=this_orbit
               if this_bpm=='BPMSW.1R5.B1':
                   orb1r5b1=this_orbit
               if this_bpm=='BPMS.2L5.B1':
                   orb2l5b1=this_orbit
               if this_bpm=='BPMS.2R5.B1':
                   orb2r5b1=this_orbit
            
               if this_bpm=='BPMSW.1L1.B1':
                   orb1l1b1=this_orbit
               if this_bpm=='BPMSW.1R1.B1':
                   orb1r1b1=this_orbit
               if this_bpm=='BPMS.2L1.B1':
                   orb2l1b1=this_orbit
               if this_bpm=='BPMS.2R1.B1':
                   orb2r1b1=this_orbit

               if this_bpm=='BPMSW.1L5.B2':
                   orb1l5b2=this_orbit
               if this_bpm=='BPMSW.1R5.B2':
                   orb1r5b2=this_orbit
               if this_bpm=='BPMS.2L5.B2':
                   orb2l5b2=this_orbit
               if this_bpm=='BPMS.2R5.B2':
                   orb2r5b2=this_orbit
            
               if this_bpm=='BPMSW.1L1.B2':
                   orb1l1b2=this_orbit
               if this_bpm=='BPMSW.1R1.B2':
                   orb1r1b2=this_orbit
               if this_bpm=='BPMS.2L1.B2':
                   orb2l1b2=this_orbit
               if this_bpm=='BPMS.2R1.B2':
                   orb2r1b2=this_orbit
               ############################################### ### if the bpm is located in the arc, and has a non-nan value of the orbit, append to a list of the orbit in the beam1/beam2 arcs.
               if not is_arc_BPM(this_bpm):
                   continue
               if math.isnan(this_orbit):
                   continue
               if re.search('\.B1',this_bpm):
                   arc_orbit_b1.append(this_orbit)
               if re.search('\.B2',this_bpm):
                   arc_orbit_b2.append(this_orbit)
           
           ################################################### ### BEAM1: check that have sufficient data to calculate the orbit in the arcs. if so calculate the mean and rms closed orbit in the arcs
           if len(arc_orbit_b1)==0:
               if b1_arc_nan==False:
                   print 'Warning: No Beam1 arc orbit data found. Arc orbit properties set to NaN as of '+convert_to_data_output_format(time)
                   b1_arc_nan=True
               mean_arc_orbit_b1=float('NaN')
               rms_arc_orbit_b1=float('NaN')
           elif len(arc_orbit_b1)>0 and len(arc_orbit_b1)<20:
               if b1_arc_nan==False:
                   print 'Warning: <20 Arc BPMs were found for Beam1. Arc orbit properties set to NaN as of '+convert_to_data_output_format(time)
                   b1_arc_nan=True
               mean_arc_orbit_b1=float('NaN')
               rms_arc_orbit_b1=float('NaN')
           else:
               if b1_arc_nan==True:
                   print 'Warning: As of '+convert_to_data_output_format(time)+' now have sufficient data to calculate Beam1 arc orbit.'
                   b1_arc_nan=False
               mean_arc_orbit_b1=np.mean(np.array(arc_orbit_b1))
               rms_arc_orbit_b1=np.sqrt(np.mean(pow(np.array(arc_orbit_b1),2)))
           #################################################### ### BEAM2: check that have sufficient data to calculate the orbit in the arcs. if so calculate the mean and rms closed orbit in the arcs
           if len(arc_orbit_b2)==0:
               if b2_arc_nan==False:
                   print 'Warning: No Beam2 arc orbit data found. Arc orbit properties set to NaN as of '+convert_to_data_output_format(time)
                   b2_arc_nan=True
               mean_arc_orbit_b2=float('NaN')
               rms_arc_orbit_b2=float('NaN')
           elif len(arc_orbit_b2)>0 and len(arc_orbit_b2)<20:
               if b2_arc_nan==False:
                   print 'Warning: <20 Arc BPMs were found for Beam2. Arc orbit properties set to NaN as of '+convert_to_data_output_format(time)
                   b2_arc_nan=True
               mean_arc_orbit_b2=float('NaN')
               rms_arc_orbit_b2=float('NaN')
           else:
               if b2_arc_nan==True:
                   print 'Warning: As of '+convert_to_data_output_format(time)+' now have sufficient data to calculate Beam1 arc orbit.'
                   b2_arc_nan=False
               mean_arc_orbit_b2=np.mean(np.array(arc_orbit_b2))
               rms_arc_orbit_b2=np.sqrt(np.mean(pow(np.array(arc_orbit_b2),2)))
           ################################################### ### try to calculate the xing angles etc. if fail set the derived values to nan
           try:
               xing_ir1b1=1e6*math.atan2((-1*orb1l1b1+orb1r1b1)*1e-3,42.95)    ### orbit values converted to meters, giving angle in rad, then converted to urad
           except:
               xing_ir1b1=float('NaN')
           try:
               xing_ir1b2=1e6*math.atan2(-1*(-1*orb1l1b2+orb1r1b2)*1e-3,42.95) ### orbit values converted to meters, giving angle in rad, then converted to urad
           except:
               xing_ir1b2=float('NaN')
           try:
               xing_ir5b1=1e6*math.atan2((-1*orb1l5b1+orb1r5b1)*1e-3,42.95)    ### orbit values converted to meters, giving angle in rad, then converted to urad  
           except:
               xing_ir5b1=float('NaN')
           try:
               xing_ir5b2=1e6*math.atan2(-1*(-1*orb1l5b2+orb1r5b2)*1e-3,42.95) ### orbit values converted to meters, giving angle in rad, then converted to urad
           except:
               xing_ir5b2=float('NaN')
           ###################################################
           orbit_properties.append([time,mean_arc_orbit_b1,mean_arc_orbit_b2,rms_arc_orbit_b1,rms_arc_orbit_b2,xing_ir1b1,xing_ir1b2,xing_ir5b1,xing_ir5b2,orb1l1b1,orb1r1b1,orb2l1b1,orb2r1b1,orb1l1b2,orb1r1b2,orb2l1b2,orb2r1b2,orb1l5b1,orb1r5b1,orb2l5b1,orb2r5b1,orb1l5b2,orb1r5b2,orb2l5b2,orb2r5b2])
       all_orbit_properties.append(orbit_properties)
    
    if len(all_orbit_properties[0])!=len(all_orbit_properties[1]):
        print 'Warning: Have different amount of orbit data for Beam1 and Beam2.'

    print '---> Merging data from the orbit variales'
    if len(all_orbit_properties)>2:
        sys.exit('Too many orbit variables. Normally expect one for H and one for V. Will have to update this part of the code to continue.')
    derived_orbit_properties,ir_bpm_orbits = merge_orbit_properties(all_orbit_properties)

    return derived_orbit_properties,ir_bpm_orbits
###########################################################################################################################################################################
def write_orbit_properties(outputpath,derived_orbit_properties,ir_bpm_orbits):

    out_derived_orbit_properties=cp.deepcopy(derived_orbit_properties) ### make local copies of the orbit lists so can convert dt objects to str for output
    out_ir_bpm_orbits=cp.deepcopy(ir_bpm_orbits)                       ### same
    
    orbit_properties_file=outputpath+'data.orbit.arc.xing.csv'
    ir_orbit_file=outputpath+'data.orbit.ir_bpms.csv'

    print '---> Writing Mean/RMS Arc orbit and crossing angles vs time to file: '+orbit_properties_file
    wout=open(orbit_properties_file,'w')
    csvwout=csv.writer(wout,delimiter=',')
    for t in range(len(out_derived_orbit_properties)):
        if t>0:
            out_derived_orbit_properties[t][0]=convert_to_data_output_format(out_derived_orbit_properties[t][0])
        csvwout.writerow(out_derived_orbit_properties[t])
    wout.close()

    print '---> Writing orbit at IR BPMs vs time to file: '+ir_orbit_file
    wout=open(ir_orbit_file,'w')
    csvwout=csv.writer(wout,delimiter=',')
    for t in range(len(out_ir_bpm_orbits)):
        if t>0:
            out_ir_bpm_orbits[t][0]=convert_to_data_output_format(out_ir_bpm_orbits[t][0])
        csvwout.writerow(out_ir_bpm_orbits[t])
    wout.close()
    return
###########################################################################################################################################################################
def read_orbit_properties(datapath):
    derived_orbit_properties=[]
    ir_bpm_orbits=[]
    orbit_properties_file=datapath+'data.orbit.arc.xing.csv'
    ir_orbit_file=datapath+'data.orbit.ir_bpms.csv'

    print '---> Reading Mean/RMS Arc orbit and crossing angles vs time from file: '+orbit_properties_file
    rdat=open(orbit_properties_file,'r')
    csvrdat=csv.reader(rdat,delimiter=',',skipinitialspace=True)
    for row in csvrdat:
        this_data=row
        if re.search('#',row[0]):
            pass
        else:
            this_data=row
            this_data[0]=convert_from_data_output_format(this_data[0])
        derived_orbit_properties.append(this_data)
    rdat.close()

    print '---> Reading orbit at IR BPMs vs time from file: '+ir_orbit_file
    rdat=open(ir_orbit_file,'r')
    csvrdat=csv.reader(rdat,delimiter=',',skipinitialspace=True)
    for row in csvrdat:
        this_data=row
        if re.search('#',row[0]):
            pass
        else:
            this_data=row
            this_data[0]=convert_from_data_output_format(this_data[0])
        ir_bpm_orbits.append(this_data)
    rdat.close()

    return derived_orbit_properties,ir_bpm_orbits
###########################################################################################################################################################################
def extract_currents(outputpath,knobs,circuits_in_knob,dtstart,dtend,ldb):
    knobs_current_data=[]
    for k in range(len(knobs)):
        print 'Extracting measured currents for knob: '+knobs[i]
        current_data=[]
        for c in range(len(circuits_in_knob[i])):
            this_circuit=circuits_in_knob[i][c]
            #print '---> Looking for circuit: '+this_circuit
            search_term='%'+this_circuit+'%I_MEAS%'
            search_result=ldb.search(search_term)
            if len(search_result)>1:
                print 'Found timber variables:'
                print search_result
                sys.exit('Circuit appears to not be uniquely defined....')
            else:
                current_var=search_result[0]
            print '---> Extracting '+current_var
            timber_st = convert_to_timber_input_format(dtstart)
            timber_et = convert_to_timber_input_format(dtend)
            this_current_data = ldb.get(current_var,timber_st,timber_et,unixtime=False)
            datatimes=this_current_data[current_var][0]
            datacurrent=this_current_data[current_var][1]
            for t in range(len(datatimes)):
                somedata=[datatimes[t]]
                for cc in range(len(circuits_in_knob[i])):
                    if cc==c:
                        somedata.append(datacurrent[t])
                    else:
                        somedata.append(None)
                current_data.append(somedata)
        current_data=sorted(current_data,key=itemgetter(0))
        
        for t in range(1,len(current_data)):
            if current_data[t][0]==current_data[t-1][0]:
                for c in range(len(current_data[t])):
                    if current_data[t][c]==None:
                        current_data[t][c]=current_data[t-1][c]

        reduced_current_data=[]
        for t in range(1,len(current_data)):
            if current_data[t][0]!=current_data[t-1][0]:
                reduced_current_data.append(current_data[t-1])
            if current_data[t][0]==current_data[t-1][0] and t==len(current_data)-1:
                reduced_current_data.append(current_data[t])

        for t in range(len(reduced_current_data)):
            for c in range(len(reduced_current_data[t])):
                if reduced_current_data[t][c]==None:
                    reduced_current_data[t][c]=float('Nan')
        
        header=[]
        ###header.append('#timestamp') ### EHM 28/04
        header.append('#time')
        for c in range(len(circuits_in_knob[i])):
            header.append(circuits_in_knob[i][c])
        reduced_current_data.insert(0,header)
        
        knobs_current_data.append(reduced_current_data)
    return knobs_current_data
###########################################################################################################################################################################
def write_currents(outputpath,knobs,knobs_current_data):
    for i in range(len(knobs)):
        knobname_for_file = convert_name(knobs[i])
        this_current_data=cp.deepcopy(knobs_current_data[i]) ### make local copy of the current data for this particular knob. use deepcopy so can convert dt objects to string output.
        ifile=outputpath+'data.Imeas.'+knobname_for_file+'.csv'
        print '---> Writing measured currents to file:  '+ifile
        wout=open(ifile,'w')
        csvwout=csv.writer(wout,delimiter=',')
        for t in range(len(this_current_data)):
            if t==0 and re.search('#',this_current_data[t][0]):
                csvwout.writerow(this_current_data[t])
            else:
                this_current_data[t][0]=convert_to_data_output_format(this_current_data[t][0])
                csvwout.writerow(this_current_data[t])
        wout.close()
    return
###########################################################################################################################################################################
def read_currents(datapath,knobs):
    knobs_current_data=[]
    for i in range(len(knobs)):
        this_current_data=[]
        knobname_for_file = convert_name(knobs[i])
        ifile=datapath+'data.Imeas.'+knobname_for_file+'.csv'
        if os.path.exists(ifile):
            print '---> Reading measured currents from file:  '+ifile
        else:
            sys.exit('Could not find existing current data: '+ifile)
        rdat=open(ifile,'r')
        csvrdat=csv.reader(rdat,delimiter=',',skipinitialspace=True)
        for row in csvrdat:
             if re.search('#',row[0]):
                 this_current_data.append(row)
             else:
                 thisdata=row
                 thisdata[0]=convert_from_data_output_format(thisdata[0])
                 this_current_data.append(thisdata)
        rdat.close()
        if len(this_current_data)==0:
            sys.exit('Found data file, but appears to be empty... Stopping.')
        knobs_current_data.append(this_current_data)
    return knobs_current_data

###########################################################################################################################################################################
###########################################################################################################################################################################
def extract_BBQ(outputpath,dtstart,dtend,ldb):
    
    bbq_vars=[]
    search_term='%EIGEN%FREQ%'
    search_result=ldb.search(search_term)
    for e in range(len(search_result)):
        if re.search('BOFSU:',search_result[e]) or re.search('BQBBQ.CONTINUOUS',search_result[e]):
            bbq_vars.append(search_result[e])
    
    search_term='%COUPL%ABS%'
    search_result=ldb.search(search_term)
    for e in range(len(search_result)):
        if re.search('BOFSU:',search_result[e]) or re.search('BQBBQ.CONTINUOUS',search_result[e]):
            bbq_vars.append(search_result[e])

    timber_st = convert_to_timber_input_format(dtstart)
    timber_et = convert_to_timber_input_format(dtend)

    temp_bbq_data=[]
    for e in range(len(bbq_vars)):
        bbq_var=bbq_vars[e]
        print '---> Extracting '+bbq_var
        this_bbq_data = ldb.get(bbq_var,timber_st,timber_et,unixtime=False)
        datatimes=this_bbq_data[bbq_var][0]
        datacurrent=this_bbq_data[bbq_var][1]
        for t in range(len(datatimes)):
            this_time=datatimes[t]
            this_value=datacurrent[t]
            this_data=[]
            this_data.append(this_time)
            for z in range(len(bbq_vars)):
                if z==e:
                    this_data.append(this_value)
                else:
                    this_data.append(None)
            temp_bbq_data.append(this_data)
    
    sorted_bbq_data=sorted(temp_bbq_data,key=itemgetter(0))
    for m in range(1,len(sorted_bbq_data)):
        if sorted_bbq_data[m][0]==sorted_bbq_data[m-1][0]:
            for n in range(len(sorted_bbq_data[m])):
                if sorted_bbq_data[m][n]==None:
                    sorted_bbq_data[m][n]=sorted_bbq_data[m-1][n]

    bbq_data=[]
    for m in range(1,len(sorted_bbq_data)):
        if sorted_bbq_data[m][0]!=sorted_bbq_data[m-1][0]:
            bbq_data.append(sorted_bbq_data[m-1])
        if sorted_bbq_data[m][0]==sorted_bbq_data[m-1][0] and m==len(sorted_bbq_data):
            bbq_data.append(sorted_bbq_data[m])

    return bbq_vars,bbq_data
###########################################################################################################################################################################
def write_bbq(outputpath,bbq_vars,bbq_data):
    this_bbq_data=cp.deepcopy(bbq_data) ### deepcopy bbq data so that when convert dt objects to timestamps don't edit the original array
    bbqfile=outputpath+'data.BBQ.csv'
    print '---> Writing BBQ data to file: '+bbqfile
    wout=open(bbqfile,'w')
    csvwout=csv.writer(wout,delimiter=',')
    
    header=['#date']
    for v in range(len(bbq_vars)):
        header.append(bbq_vars[v])
    
    csvwout.writerow(header)
    for t in range(len(this_bbq_data)):
        this_bbq_data[t][0]=convert_to_data_output_format(this_bbq_data[t][0])
        for n in range(len(this_bbq_data[t])):
            if this_bbq_data[t][n]==None:
                this_bbq_data[t][n]='nan'
        csvwout.writerow(this_bbq_data[t])
    wout.close()
    return

###########################################################################################################################################################################
def write_currents(outputpath,knobs,knobs_current_data):
    for i in range(len(knobs)):
        knobname_for_file = convert_name(knobs[i])
        this_current_data=cp.deepcopy(knobs_current_data[i]) ### make local copy of the current data for this particular knob. use deepcopy so can convert dt objects to string output.
        ifile=outputpath+'data.Imeas.'+knobname_for_file+'.csv'
        print '---> Writing measured currents to file:  '+ifile
        wout=open(ifile,'w')
        csvwout=csv.writer(wout,delimiter=',')
        for t in range(len(this_current_data)):
            if t==0 and re.search('#',this_current_data[t][0]):
                csvwout.writerow(this_current_data[t])
            else:
                this_current_data[t][0]=convert_to_data_output_format(this_current_data[t][0])
                csvwout.writerow(this_current_data[t])
        wout.close()
    return
#######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
def main_get_data(extract,knobs,circuits_in_knob,all_data_platteaus,dtstart,dtend):
    ldb = pytimber.LoggingDB()
    print '\n\n'
    print 'Obtaining BPM orbit data'
    print '-------------------------'
    orbitvars=define_orbit_vars()
    if extract==True:
        all_meta_to_use,all_orbit_data = extract_orbit(outputpath,orbitvars,dtstart,dtend,ldb)
        
    else:
        all_meta_to_use,all_orbit_data = read_orbit(datapath,orbitvars)

    if extract==True:
        derived_orbit_properties,ir_bpm_orbits = get_inferred_orbit(orbitvars,all_meta_to_use,all_orbit_data)
        write_orbit_properties(outputpath,derived_orbit_properties,ir_bpm_orbits)
    else:
        derived_orbit_properties,ir_bpm_orbits = read_orbit_properties(datapath)
    
    print '\n\n'
    print 'Obtaining Measured currents'
    print '----------------------------'
    if extract==True:    
        knobs_current_data = extract_currents(outputpath,knobs,circuits_in_knob,dtstart,dtend,ldb)
        write_currents(outputpath,knobs,knobs_current_data)
    else:
        knobs_current_data = read_currents(datapath,knobs)


    print '\n\n'
    print 'Obtaining BBQ data'
    print '-------------------'
    bbq_vars,bbq_data = extract_BBQ(outputpath,dtstart,dtend,ldb)
    write_bbq(outputpath,bbq_vars,bbq_data)

    #search_terms='%BOFSU%position%'
    #search_results=ldb.search(search_terms)
    #print search_results
    return
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################
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
    knobs,dtstart,dtend,outputpath,extract,datapath,acdipole,OFB = parse_options()
    ####################################################################
    print '\n\n'
    print '----------------------------------------------------------------------------------'

    print 'Will analyse orbit bump scan between '+dtstart.strftime('%Y-%m-%d %H:%M:%S.%f')+' and '+dtend.strftime('%Y-%m-%d %H:%M:%S.%f')+''
    for i in range(len(knobs)):
        print 'Will analyse knob: '+knobs[i]
    print ''
    if extract==True:
        print 'Will extract data from timber'
    else:
        print 'Will use existing data located at: '+datapath
    print ''
    if acdipole==True:
        print 'Will search for AC-dipole excitations during the scan'
    else:
        print 'Will NOT search for AC-dipole excitations during the scan'
    print ''
    print 'Output will be written to: '+outputpath
    print '----------------------------------------------------------------------------------\n\n'
    ####################################################################
    circuits_in_knob,all_data_platteaus = main_get_platteaus(knobs,dtstart,dtend,extract,outputpath,datapath,OFB)
    ####################################################################
    main_get_data(extract,knobs,circuits_in_knob,all_data_platteaus,dtstart,dtend)
    ####################################################################
    main_analysis_bbq()
    ####################################################################
    main_analysis_acd()
    ####################################################################
    main_model()
    ####################################################################
    ####################################################################
    ####################################################################

