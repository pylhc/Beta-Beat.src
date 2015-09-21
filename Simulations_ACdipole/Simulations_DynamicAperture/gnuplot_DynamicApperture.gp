set terminal postscript enhanced color solid 15 font ",20"
set output "./plot_dynamic_apperture.eps"

set notitle
set xlabel "{/Symbol A}_{x} [mm]
set ylabel "{/Symbol A}_{y} [mm]
set xtics 0.4
set ytics 0.4
set xrange [0.3:3.7]
set yrange [0.3:3.7]
set view map
set key 3.4,4.4
set key invert
set key spacing 1.5

set size 1.4,1
set origin 0,0
set multiplot

set size 0.666,1
set origin 0,0

set logscale cb
set cbrange [1e-4:5]

splot "dynamic_apperture_single.dat" u 1:2:($5!=0?$5*1e3:1/0) with points pointtype 5 pointsize 3 lt rgb "grey" linewidth 30 t "single kick: stable", "DynamicApperture/dynamic_apperture.dat" u 1:2:($5!=0?$5*1e3:1/0) with points pointtype 5 pointsize 3 palette linewidth 30 t "AC dipole: P2P_{x} after ramp down [mm]"

set size 0.666,1
set origin 0.7,0

set logscale cb
set cbrange [1e-4:1.5]

splot "dynamic_apperture_single.dat" u 1:2:($6!=0?$6*1e3:1/0) with points pointtype 5 pointsize 3 lt rgb "grey" linewidth 30 t "single kick: stable", "DynamicApperture/dynamic_apperture.dat" u 1:2:($6!=0?$6*1e3:1/0) with points pointtype 5 pointsize 3 palette linewidth 30 t "AC dipole: P2P_{x} after ramp down [mm]"

unset multiplot


