#!/bin/bash

gnuplot << TOEND

set term postscript eps color enhanced "Helvetica" 20

set output 'light.eps'
set grid
set title "Light field, {/Symbol l}"
set ylabel "{/Symbol l}
set xlabel "t"
set yrange [0:250]

plot \
'data.txt' using 1:2 title '|{/Symbol l}|^2' w lines

set output 'phi.eps'
set grid
set title "Momentum Matrix, {/Symbol F}"
set ylabel "{/Symbol F}
set xlabel "t"
set yrange [-1:0.5]

plot \
'data.txt' using 1:3 title 'S_{Z}' w lines
TOEND

convert -density 150 light.eps light.png
convert -density 150 phi.eps phi.png
rm light.eps
rm phi.eps
