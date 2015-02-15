#!/bin/bash
cat << TOEND > plotfile.gp

set term postscript eps color enhanced "TimesNewRoman" 28 dashlen 4

set output 'OmTil_2.500-NSweep.eps'
set label '~{/Symbol w}{0.6\~} = 2.1' at graph 0.8, 0.9
set ylabel "|{/Symbol l}|^2, cavity field"
set xlabel "P, pump field"
set yrange [0:3]
set xrange [0.4:0.65]
set xtics 0.1
set format x "%.2f"
set key top left

a = 'OmTil_2.500-n_1.dat'
b = 'OmTil_2.500-n_2.dat'
c = 'OmTil_2.500-n_3.dat'
e = 'OmTil_2.500-n_5.dat'
d = 'OmTil_2.500-n_10.dat'

set palette defined (0 "blue",1 "red", 2 "yellow")
set pal mod HSV
set pal def (0 0.2 1 1, 1 1 1 1)

set cbrange [0 to 10]
unset colorbox

plot \
a w l title "n = 1" lw 7 lc palette cb 1, \
b w l title "n = 2" lw 7 lc palette cb 2 , \
c w l title "n = 3" lw 7 lc palette cb 3, \
e w l title "n = 5" lw 7 lt 5 lc palette cb 5, \
d w l title "n = 10" lw 7 lt 4 lc palette cb 10
TOEND

gnuplot plotfile.gp

convert -density 300 OmTil_2.500-NSweep.eps OmTil_2.500-NSweep.png 
