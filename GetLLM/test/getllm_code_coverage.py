'''
Created on 10 Apr 2013

@author: vimaier

@version: 1.0.0

This module runs several times GetLLM  with different input data to get the complete code coverage.

To get the code coverage use the tool 'Code coverage' from PyDev:
http://pydev.org/manual_adv_coverage.html

Change history:
 - 
'''


import GetLLM

# Global variables
my_outputpath = ""
my_files_to_analyse = ""
my_twiss_model_file = ""
    
my_accel = "LHCB1"
my_dict_file = "0"
my_lhcphase = "0"
my_BPMU = "um"
my_COcut = "4000"
my_NBcpl = 2
my_TBTana = "SUSSIX"
my_higher_order ="1"

def run_getllm():
    GetLLM.main(
                outputpath=my_outputpath,
                dict_file=my_dict_file,
                files_to_analyse=my_files_to_analyse,
                twiss_model_file=my_twiss_model_file,
                accel=my_accel,
                lhcphase=my_lhcphase,
                BPMU=my_BPMU,
                COcut=float(my_COcut),
                NBcpl=my_NBcpl,
                TBTana=my_TBTana,
                higher_order=my_higher_order) 

if __name__ == "__main__":
    
    
    #getsuper_all_three
    print "========================================================================================"
    print "getsuper_all_three"
    print "========================================================================================"
    my_outputpath = "data/getsuper_all_three/output/to_check"
    my_files_to_analyse = "data/getsuper_all_three/input/src_files/bpms1.gz,"
    my_files_to_analyse += "data/getsuper_all_three/input/src_files/bpms2.gz,"
    my_files_to_analyse += "data/getsuper_all_three/input/src_files/bpms3.gz"
    my_twiss_model_file = "data/getsuper_all_three/input/model/twiss.dat"
    run_getllm()  
      
    #getsuper_all_three with nbcpl = 1
    print "========================================================================================"
    print "getsuper_all_three with nbcpl = 1"
    print "========================================================================================"
    my_BPMU = "um"
    
    my_NBcpl = 1
    run_getllm()   
    my_NBcpl = 2
          
    #getsuper_data1
    print "========================================================================================"
    print "getsuper_data1"
    print "========================================================================================"
    my_outputpath = "data/getsuper_data1/output/to_check"
    my_files_to_analyse = "data/getsuper_data1/input/src_files/bpms.gz"
    my_twiss_model_file = "data/getsuper_data1/input/model/twiss.dat"
    run_getllm()        
          
    #getsuper_data2
    print "========================================================================================"
    print "getsuper_data2"
    print "========================================================================================"
    my_outputpath = "data/getsuper_data2/output/to_check"
    my_files_to_analyse = "data/getsuper_data2/input/src_files/bpms.gz"
    my_twiss_model_file = "data/getsuper_data2/input/model/twiss.dat"
    run_getllm()        
          
    #getsuper_data3
    print "========================================================================================"
    print "getsuper_data3"
    print "========================================================================================"
    my_outputpath = "data/getsuper_data3/output/to_check"
    my_files_to_analyse = "data/getsuper_data3/input/src_files/bpms.gz"
    my_twiss_model_file = "data/getsuper_data3/input/model/twiss.dat"
    run_getllm()        

    #run_from_gui with my_BPMU = "mm"
    print "========================================================================================"
    print "run_from_gui with my_BPMU = \"mm\""
    print "========================================================================================"
    my_outputpath = "data/run_from_gui/output/to_check"
    my_files_to_analyse = "data/run_from_gui/input/src_files/ALLBPMs"
    my_twiss_model_file = "data/run_from_gui/input/model/twiss.dat"
    my_BPMU = "mm"
    run_getllm()        
    my_BPMU = "um"
          
    #run_with_ac
    print "========================================================================================"
    print "run_with_ac"
    print "========================================================================================"
    my_outputpath = "data/run_with_ac/output/to_check"
    my_files_to_analyse = "data/run_with_ac/input/src_files/ALLBPMs"
    my_twiss_model_file = "data/run_with_ac/input/model/twiss.dat"
    run_getllm()        
          
    #run_with_ip
    print "========================================================================================"
    print "run_with_ip"
    print "========================================================================================"
    my_outputpath = "data/run_with_ip/output/to_check"
    my_files_to_analyse = "data/run_with_ip/input/src_files/ALLBPMs"
    my_twiss_model_file = "data/run_with_ip/input/model/twiss.dat"
    run_getllm()        
          
    #run_with_ip_and_ac
    print "========================================================================================"
    print "run_with_ip_and_ac"
    print "========================================================================================"
    my_outputpath = "data/run_with_ip_and_ac/output/to_check"
    my_files_to_analyse = "data/run_with_ip_and_ac/input/src_files/ALLBPMs"
    my_twiss_model_file = "data/run_with_ip_and_ac/input/model/twiss.dat"
    run_getllm()        
         
    #run1
    print "========================================================================================"
    print "run1"
    print "========================================================================================"
    my_outputpath = "data/run1/output/to_check"
    my_files_to_analyse = "data/run1/input/src_files/ALLBPMs"
    my_twiss_model_file = "data/run1/input/model/twiss.dat"
    run_getllm()        
    
     

          
              
