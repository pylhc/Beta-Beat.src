mkdir DynamicApperture

for amp_x in 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6
do
    for amp_y in 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6
    do

	mkdir DynamicApperture/Amplitude_$amp_x-$amp_y
	cp job.Beam1.ACdipoletracking.6500GEV.40cm.madx DynamicApperture/Amplitude_$amp_x-$amp_y/job.$amp_x-$amp_y.mask 
	sed -i "s/XAMPLITUDE/$amp_x/g" DynamicApperture/Amplitude_$amp_x-$amp_y/job.$amp_x-$amp_y.mask 
	sed -i "s/YAMPLITUDE/$amp_y/g" DynamicApperture/Amplitude_$amp_x-$amp_y/job.$amp_x-$amp_y.mask 
	bsub -o ./directory_temp/ -q 8nh "/afs/cern.ch/eng/lhc_online_model/pro/virtualcorrectors/madx64 < PATH_TO_FOLDER/DynamicApperture/Amplitude_$amp_x-$amp_y/job.$amp_x-$amp_y.mask"

    done
done


