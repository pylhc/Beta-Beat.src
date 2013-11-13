'''
Created on ???

@author: ???

@version: 1.0.1

getdiff.py will be executed after GetLLM.

Provide as first argument the path to the output files of GetLLM.

getdiff needs as input:
    - getbetax_free.out or getbetax.out
    - getbetay_free.out or getbetay.out
    - getcouple_free.out or getcouple.out
    - getDx.out
    - getDy.out
    - twiss_cor.dat
    - twiss_no.dat

getdiff produces the following output in the directory stated as argument:
    - bbx.out
    - bby.out
    - couple.out
    - dx.out
    - dy.out


Change history:
 - 1.0.1, vimaier, 04th June 2013:
     Added module description
     Reformatted into sections parse_args(), main() and main invocation
     Renamed variables
     Removed warnings from static analysis
'''


import os
import sys

import __init__  # @UnusedImport __init__ adds the path to Beta-Beat.src
import Python_Classes4MAD.metaclass





#===================================================================================================
# parse_args()-function
#===================================================================================================
def parse_args():
    ''' Parses the sys.argv[1], checks for valid input and returns the path to src files  '''
    try:
        path_to_src_files = sys.argv[1]
    except IndexError:
        print >> sys.stderr, str("Provide path to GetLLM output as first command line argument. E.g:"+
                                    "  python getdiff.py /path/to/getllm/output/")
        sys.exit(1)

    if not os.path.isdir(path_to_src_files):
        print >> sys.stderr, "No valid directory:", path_to_src_files
    return path_to_src_files

#===================================================================================================
# main()-function
#===================================================================================================
def main(path):
    '''
    :Parameters:
        'path': string
            Path to src files. Will also be used for output files.
    :Return: int
        0 if execution was successful otherwise !=0
    '''
    twiss_cor = Python_Classes4MAD.metaclass.twiss(os.path.join(path,'twiss_cor.dat'))
    twiss_no = Python_Classes4MAD.metaclass.twiss(os.path.join(path,'twiss_no.dat'))
    twiss_cor.Cmatrix()
    twiss_no.Cmatrix()


    # normal quad

    file_bbx = open( os.path.join(path,"bbx.out"), "w")
    file_bby = open( os.path.join(path,"bby.out"), "w")

    print >> file_bbx,"NAME S MEA ERROR MODEL"
    print >> file_bbx,"%s %le %le %le %le"

    print >> file_bby,"NAME S MEA ERROR MODEL"
    print >> file_bby,"%s %le %le %le %le"

    if os.path.exists( os.path.join(path,'getbetax_free.out') ):
        twiss_getbetax = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getbetax_free.out') )
    else:
        twiss_getbetax = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getbetax.out') )

    for i in range(len(twiss_getbetax.NAME)):
        bpm_name = twiss_getbetax.NAME[i]
        bpm_included = True
        try:
            check = twiss_cor.NAME[twiss_cor.indx[bpm_name]]  # @UnusedVariable
        except:
            print "No ", bpm_name
            bpm_included = False
        if bpm_included:
            j=twiss_cor.indx[bpm_name]
            t_x = twiss_getbetax # Variable for abbreviation
            print>> file_bbx,bpm_name, t_x.S[i], (t_x.BETX[i]-t_x.BETXMDL[i])/t_x.BETXMDL[i], t_x.STDBETX[i]/t_x.BETXMDL[i],(twiss_cor.BETX[j]-twiss_no.BETX[j])/twiss_no.BETX[j]

    if os.path.exists( os.path.join(path,'getbetay_free.out') ):
        twiss_getbetay = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getbetay_free.out') )
    else:
        twiss_getbetay = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getbetay.out') )

    for i in range(len(twiss_getbetay.NAME)):
        bpm_name=twiss_getbetay.NAME[i]
        bpm_included = True
        try:
            check = twiss_cor.NAME[twiss_cor.indx[bpm_name]] # @UnusedVariable
        except:
            print "No ", bpm_name
            bpm_included = False
        if bpm_included:
            j=twiss_cor.indx[bpm_name]
            t_y = twiss_getbetay # Variable for abbreviation
            print>> file_bby,bpm_name, t_y.S[i], (t_y.BETY[i]-t_y.BETYMDL[i])/t_y.BETYMDL[i], t_y.STDBETY[i]/t_y.BETYMDL[i],(twiss_cor.BETY[j]-twiss_no.BETY[j])/twiss_no.BETY[j]

    file_bbx.close()
    file_bby.close()

    try:
        twiss_getdx = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getDx.out') )

        file_dx = open( os.path.join(path,"dx.out") ,"w")

        print >> file_dx,"NAME S MEA ERROR MODEL"
        print >> file_dx,"%s %le %le %le %le"

        for i in range(len(twiss_getdx.NAME)):
            bpm_name = twiss_getdx.NAME[i]
            bpm_included = True
            try:
                check = twiss_cor.NAME[twiss_cor.indx[bpm_name]] # @UnusedVariable
            except:
                print "No ", bpm_name
                bpm_included = False
            if bpm_included:
                j=twiss_cor.indx[bpm_name]
                print>> file_dx,bpm_name, twiss_getdx.S[i], (twiss_getdx.DX[i]-twiss_getdx.DXMDL[i]), twiss_getdx.STDDX[i],(twiss_cor.DX[j]-twiss_no.DX[j])

        file_dx.close()
    except IOError:
        print "NO dispersion"
    except AttributeError:
        print "Empty table in getDx.out?! NO dispersion"


    # skew quad

    file_couple = open( os.path.join(path,"couple.out") ,"w")


    print >> file_couple,"NAME S F1001re F1001im F1001e F1001re_m F1001im_m"
    print >> file_couple,"%s %le %le %le %le %le %le"



    if os.path.exists( os.path.join(path,'getcouple_free.out') ):
        twiss_getcouple = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getcouple_free.out') )
    else:
        twiss_getcouple = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getcouple.out') )


    for i in range(len(twiss_getcouple.NAME)):
        bpm_name = twiss_getcouple.NAME[i]
        bpm_included = True
        try:
            check = twiss_cor.NAME[twiss_cor.indx[bpm_name]] # @UnusedVariable
        except:
            print "No ", bpm_name
            bpm_included = False
        if bpm_included:
            j = twiss_cor.indx[bpm_name]
            print >> file_couple,bpm_name, twiss_getcouple.S[i],twiss_getcouple.F1001R[i],twiss_getcouple.F1001I[i],twiss_getcouple.FWSTD1[i],twiss_cor.f1001[j].real,twiss_cor.f1001[j].imag




    try:
        twiss_getdy = Python_Classes4MAD.metaclass.twiss( os.path.join(path,'getDy.out') )

        file_dy = open( os.path.join(path,"dy.out") ,"w")

        print >> file_dy,"NAME S MEA ERROR MODEL"
        print >> file_dy,"%s %le %le %le %le"

        for i in range(len(twiss_getdy.NAME)):
            bpm_name = twiss_getdy.NAME[i]
            bpm_included = True
            try:
                check = twiss_cor.NAME[twiss_cor.indx[bpm_name]] # @UnusedVariable
            except:
                print "No ", bpm_name
                bpm_included = False
            if bpm_included:
                j = twiss_cor.indx[bpm_name]
                print>> file_dy,bpm_name, twiss_getdy.S[i], (twiss_getdy.DY[i]-twiss_getdy.DYMDL[i]), twiss_getdy.STDDY[i],(twiss_cor.DY[j]-twiss_no.DY[j])

        file_dy.close()
    except IOError:
        print "NO dispersion."
    except AttributeError:
        print "Empty table in getDy.out?! NO dispersion"

    return 0


#===================================================================================================
# main invocation
#===================================================================================================
def _start():
    path_to_src_files = parse_args()

    print "Start getdiff.main..."
    return_value = main(path = path_to_src_files)
    print "getdiff.main finished with",return_value

    sys.exit(return_value)

if __name__ == "__main__":
    _start()
