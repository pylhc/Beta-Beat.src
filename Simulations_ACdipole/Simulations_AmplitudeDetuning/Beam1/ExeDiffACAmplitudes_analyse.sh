constAmp=0.5

for amp in 0.5 0.75 1.0 1.25 1.5
do
    cp AddBpmErrors.sh DetuningWithAmplitude_horizontal/Amplitude_$amp
    cp DriveFiles/* DetuningWithAmplitude_horizontal/Amplitude_$amp

    cd DetuningWithAmplitude_horizontal/Amplitude_$amp
    ./AddBpmErrors.sh
#    rm trackone
    ./command_Drive
    python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/GetLLM/GetLLM.py --accel=LHCB1 --model=twiss_elements.dat --files="ALLBPMs"  --output="./" --tbtana=SUSSIX --bpmu=m --lhcphase=1
    cd ../../
done


for amp in 0.5 0.75 1.0 1.25 1.5
do
    cp AddBpmErrors.sh DetuningWithAmplitude_vertical/Amplitude_$amp
    cp DriveFiles/* DetuningWithAmplitude_vertical/Amplitude_$amp

    cd DetuningWithAmplitude_vertical/Amplitude_$amp
    ./AddBpmErrors.sh
#    rm trackone
    ./command_Drive
    python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/GetLLM/GetLLM.py --accel=LHCB1 --model=twiss_elements.dat --files="ALLBPMs"  --output="./" --tbtana=SUSSIX --bpmu=m --lhcphase=1
    cd ../../
done

cp DataAnalysis/python_getacTuneAndAction.py DetuningWithAmplitude_horizontal
cd DetuningWithAmplitude_horizontal
python2.6 python_getacTuneAndAction.py
cd ../

cp DataAnalysis/python_getacTuneAndAction.py DetuningWithAmplitude_vertical
cd DetuningWithAmplitude_vertical
python2.6 python_getacTuneAndAction.py
cd ../

cd DataAnalysis
gnuplot gnuplot_detWithAmplitude.gp
cd ../