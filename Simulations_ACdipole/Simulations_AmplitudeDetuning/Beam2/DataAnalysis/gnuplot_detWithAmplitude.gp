set terminal postscript enhanced color solid 10 font ",12"
set encoding iso_8859_1
set key box

set output "plot.detWithAmplitude_tuneFromGetLLM_action_ac.eps"
set size 1,1
set origin 0,0
set multiplot

set yrange [0.309:0.3105]

set key top left
set size 0.5,0.5
set origin 0,0.5
set xlabel "2J_{x} [nm]"
set xrange [0:9]
set ylabel "Q_{x}"
plot "../DetuningWithAmplitude_horizontal/det_with_amp_ac.dat" u ($5*1e9):($1):($6*1e9):2 w xyerrorbars t "model" lt rgb "blue", "../DetuningWithAmplitude_horizontal/det_with_amp_ac.dat" u ($5*1e9):($1) w l t "" lt rgb "blue"

set nokey

set size 0.5,0.5
set origin 0.5,0.5
set xlabel "2J_{y} [nm]"
set xrange [0:7]
set ylabel "Q_{x}"
plot "../DetuningWithAmplitude_vertical/det_with_amp_ac.dat" u (($7*1e9<8.5)?$7*1e9:1/0):($1):($8*1e9):2 w xyerrorbars t "model" lt rgb "blue", "../DetuningWithAmplitude_vertical/det_with_amp_ac.dat" u (($7*1e9<8.5)?$7*1e9:1/0):($1) w l t "" lt rgb "blue"

set yrange [0.3185:0.3205]

set size 0.5,0.5
set origin 0,0
set xlabel "2J_{x} [nm]"
set xrange [0:9]
set ylabel "Q_{y}"
plot "../DetuningWithAmplitude_horizontal/det_with_amp_ac.dat" u ($5*1e9):($3):($6*1e9):4 w xyerrorbars t "model" lt rgb "blue", "../DetuningWithAmplitude_horizontal/det_with_amp_ac.dat" u ($5*1e9):($3) w l t "" lt rgb "blue"

set size 0.5,0.5
set origin 0.5,0
set xlabel "2J_{y} [nm]"
set xrange [0:7]
set ylabel "Q_{y}"
plot "../DetuningWithAmplitude_vertical/det_with_amp_ac.dat" u (($7*1e9<8.5)?$7*1e9:1/0):($3):($8*1e9):4 w xyerrorbars t "model" lt rgb "blue", "../DetuningWithAmplitude_vertical/det_with_amp_ac.dat" u (($7*1e9<8.5)?$7*1e9:1/0):($3) w l t "" lt rgb "blue"


unset multiplot

