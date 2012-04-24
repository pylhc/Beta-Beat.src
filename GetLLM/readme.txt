GetLLM = Get Linear Lattice function and More

Originally created by Masa. Aiba
11/Feb/2008


## Python script to obtain Linear Lattice function and More -> GetLLM
## Version 1.51
## Version-up history:V1.0, 11/Feb/2008 by Masa. Aiba
##                    V1.1, 18/Feb/2008 Debugging, add model phase and tunes to output
##                                      add function to obtain DY
##                                      add chromatic parameter (phase for non zero DPP)
##                    V1.2, 22/Feb/2008 test version for beta with all BPM
##                    V1.3, 29/Feb/2008 beta from phases is improved, averaging beta1, 2 and 3
##                    V1.31, 12/Mar/2008 debugged on alpha3
##                    V1.4, 12/Mar/2008 modify output to fit latest TSF format and to meet requests from Rogelio
##                                      fix buggs in r.m.s. beta-beat and added to the output of getbetax/y.out
##                    V1.5 Update to option parser, include BPMdictionary to filter BPMs not in Model
##                                  Rogelio, 13 March 2008
##                    V1.51, 13/Mar/2008 Modify output to fit latest TSF format again. Add STD to beta.
##                    V1.6, 15/Jul/2008 Add the integer part of tunes - assuming that the phase advance is always less than 1.0.
##                    V1.71 27/Jul/2008 Add GetCO. Filter in dispersion calculation to exclude bad bpms.
##                    V1.8, 13/Aug/2008 Add GetCoupling. - Ref. note by A. Franchi, R. T. Garcia, G. Vanbavinckhove
##                                      "Computation of the Coupling Resonance Driving term f1001 and the coupling coefficient C
##                                       from turn-by-turn single-BPM data", 28/May/2008
##                                      The GetCoupling.py is initiated by Glenn V.
##                                      and imported into GetLLM / finalized by Masa. Aiba
##                                      Some bugs are fixed - you can run even without liny file and
##                                      find the results for horizontal plane only.
##                    V1.81, 1/Sep/2008 For an accelerator in which the beam goes opposite direction to the model as in LHCB2,
##                                       the beam direction parameter (bd) is added to GetPhases.
##                                       Bug in the phi13 for the last monitor and the last but one is fixed.
##                    V1.9, 21/Oct/2008 Add the beta from spectrum height.
##                    V1.91, 08/Dec/2008 Add option - SUSSIX or SVD for input file

##                    V1.99, 13/Feb/09 Add DPX! The MEMORIAL version for Version 1.** since all the linear lattice parameters become available!
##                                     Add the option for the Harmonic analysis.
##                                     Refine coding, add several lines of comment

##                    V2.0, 17/Feb/2009 Major version up to output "More", that is, other than the linear lattice parameters.
##                                      Add off-momentum lattice (dbeta/beta)/(dp/p)
##                                      Add SPS coupling with "Pseudo-double plane BPM"-it need at this moment a list containing
##                                      pre-paired H-V monitor. This should be replaced by a clever algorithm to find good pairs automatically.
##                    V2.01, 10/Mar/2009 Fix bug on SPS double plane BPM monitor, in which missing BPM could cause an error.
##                                      Modify BetaFromAmplitude to output invariant J (by Rogelio and finalized by MA)

## Usage1 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./
## Usage2 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./


## Some rules for variable name: Dictionary is used to contain the output of function
##                               Valiable containing 'm' is a value directly obtained from measurment data
##                               Valiable containing 'mdl' is a value related to model

