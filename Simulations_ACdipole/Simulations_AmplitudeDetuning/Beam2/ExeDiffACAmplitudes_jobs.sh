constAmp=0.5

mkdir DetuningWithAmplitude_horizontal
for amp in 0.5 0.75 1.0 1.25 1.5
do
    mkdir DetuningWithAmplitude_horizontal/Amplitude_$amp
    cp job.Beam2.ACdipoletracking.6500GEV.40cm.madx DetuningWithAmplitude_horizontal/Amplitude_$amp/job.$amp.mask 
    sed -i "s/XAMPLITUDE/$amp/g" DetuningWithAmplitude_horizontal/Amplitude_$amp/job.$amp.mask 
    sed -i "s/YAMPLITUDE/$constAmp/g" DetuningWithAmplitude_horizontal/Amplitude_$amp/job.$amp.mask 
    sed -i "s/AMPLITUDE/$amp/g" DetuningWithAmplitude_horizontal/Amplitude_$amp/job.$amp.mask 
    sed -i "s/PLANE/horizontal/g" DetuningWithAmplitude_horizontal/Amplitude_$amp/job.$amp.mask 
    bsub -o ./directory_temp/ -q 8nh "/afs/cern.ch/eng/lhc_online_model/pro/virtualcorrectors/madx64 < PATH_TO_FOLDER/DetuningWithAmplitude_horizontal/Amplitude_$amp/job.$amp.mask"
done

mkdir DetuningWithAmplitude_vertical
for amp in 0.5 0.75 1.0 1.25 1.5
do
    mkdir DetuningWithAmplitude_vertical/Amplitude_$amp
    cp job.Beam2.ACdipoletracking.6500GEV.40cm.madx DetuningWithAmplitude_vertical/Amplitude_$amp/job.$amp.mask 
    sed -i "s/XAMPLITUDE/$constAmp/g" DetuningWithAmplitude_vertical/Amplitude_$amp/job.$amp.mask 
    sed -i "s/YAMPLITUDE/$amp/g" DetuningWithAmplitude_vertical/Amplitude_$amp/job.$amp.mask 
    sed -i "s/AMPLITUDE/$amp/g" DetuningWithAmplitude_vertical/Amplitude_$amp/job.$amp.mask 
    sed -i "s/PLANE/vertical/g" DetuningWithAmplitude_vertical/Amplitude_$amp/job.$amp.mask 
    bsub -o ./directory_temp/ -q 8nh "/afs/cern.ch/eng/lhc_online_model/pro/virtualcorrectors/madx64 < PATH_TO_FOLDER/DetuningWithAmplitude_vertical/Amplitude_$amp/job.$amp.mask"
done
