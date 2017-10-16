resultDir="temp"
modelDir="tests/standalone_test/models/40cmACdipole/"
measDir="tests/standalone_test/data/ADT_data/" 
python GetLLM/GetLLM.py --accel=LHCB1 --model=${modelDir}/twiss.dat --files=${measDir}/Beam1@Turn@2017_08_18@09_59_41_414.sdds.clean.new --output=$resultDir --tbtana=SUSSIX --bpmu=mm --lhcphase=1 --errordefs=$modelDir /error_deff.txt
python Correction/correct_coupleDy.py --accel=LHCB1 --path=$resultDir --cut=0.01 --errorcut=0.02,0.02 --modelcut=0.0,0.01 --rpath=/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/ --MinStr=0.000001 --Dy=1,1,0,0,0 --opt=$modelDir --Variables=coupling_knobs
python model/creator.py coupling_correction --beam=1 --nattuney=0.318 --nattunex=0.31 --accel=lhc --lhcmode=lhc_runII_2017 --output=$resultDir --optics=${modelDir}/modifiers.madx --writeto=${resultDir}/job.twiss.madx
python MODEL/LHCB/model/Corrections/getdiff.py $resultDir
python GetLLM/getllm_precision_check.py $EXCITATION --optics=40cm --analyse=sussix


