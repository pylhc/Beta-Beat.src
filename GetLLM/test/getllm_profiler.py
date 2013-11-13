'''
Created on 29 Apr 2013

@author: vimaier

@version: 1.0.0

This module uses cProfile to profile GetLLM. It produces binary files which can be read with
pstats.Stats:
http://docs.python.org/release/2.6.6/library/profile.html?highlight=pstats.stats#pstats.Stats

Change history:
 -
'''

import cProfile
import pstats
import time

import GetLLM.GetLLM  # @UnusedImport used in command string

# Global variables
my_outputpath = "data/getsuper_all_three/output/to_check"
my_files_to_analyse = "data/getsuper_all_three/input/src_files/bpms1.gz,"
my_files_to_analyse += "data/getsuper_all_three/input/src_files/bpms2.gz,"
my_files_to_analyse += "data/getsuper_all_three/input/src_files/bpms3.gz"
my_twiss_model_file = "data/getsuper_all_three/input/model/twiss.dat"

my_accel = "LHCB1"
my_dict_file = "0"
my_lhcphase = "0"
my_BPMU = "um"
my_COcut = "4000"
my_NBcpl = 2
my_TBTana = "SUSSIX"
my_higher_order ="1"

if __name__ == "__main__":

    output_name = "profiler_data/binaries/profiler_output_"+time.strftime("%Y_%m_%d")

    command = str("GetLLM.GetLLM.main(outputpath=my_outputpath,dict_file=my_dict_file,files_to_analyse=my_files_to_analyse,model_filename=my_twiss_model_file,accel=my_accel,lhcphase=my_lhcphase,BPMU=my_BPMU,COcut=float(my_COcut),NBcpl=my_NBcpl,TBTana=my_TBTana,higher_order=my_higher_order)" )
    cProfile.run(command,output_name)


    p = pstats.Stats(output_name)
    # Prints top ten internal time
    p.strip_dirs().sort_stats('time').print_stats(10)
    # Prints top ten calls count
    p.strip_dirs().sort_stats('calls').print_stats(10)
    # Prints top ten cumulative time
    p.strip_dirs().sort_stats('cumulative').print_stats(10)
    # Prints all functions sorted by name
    p.strip_dirs().sort_stats('name').print_stats()


