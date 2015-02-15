#!/bin/bash

cd fortranfiles/project/

make && ./main

gnuplot << TOEND

set term postscript eps color enhanced "Helvetica" 20

set output 'light.eps'
set grid
set title "Light field, {/Symbol l}"
set ylabel "{/Symbol l}
set xlabel "t"
set yrange [0:1]

plot \
'data1.txt' using 1:2 title 'Re({/Symbol l})' w lines, \
'data1.txt' using 1:3 title 'Im({/Symbol l})' w lines, \
'data1.txt' using 1:4 title 'abs({/Symbol l})' w lines

set output 'phi_00.eps'
set grid
set title ", {/Symbol F}"
set ylabel "{/Symbol F}
set xlabel "t"
set yrange [-5:5]

plot \
'data2.txt' using 1:2 title 'abs({/Symbol F}_00)' w lines, \
'data2.txt' using 1:3 title 'Im({/Symbol F}_00)' w lines, \
'data2.txt' using 1:4 title 'abs({/Symbol F}_00)' w lines
TOEND

convert -density 150 light.eps light.png
convert -density 150 phi_00.eps phi_00.png
rm light.eps
