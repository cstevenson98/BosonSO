#!./get.sh

dirname="Esslinger_Chaos_Omtil=1.9_P=0.469_n=10"
rt=dirname."/SI-Pump_ 0.469OmegaTil_ 1.900.dat"

dirname="Esslinger_Chaos_Omtil=0.5_P=0.261_n=10"
rt=dirname."/SI-Pump_ 0.261OmegaTil_ 0.500.dat"

ft=dirname."/power-spectrum.dat"

fname=dirname.".eps"

labstr='~{/Symbol w}{0.7\~}=1.9, {/Symbol h}=0.469, n=10'

set term postscript enhanced color eps "Times-Roman,24"
set out fname

xlab="Time, t (ms)"

unset key

set xrange [0 to 1.6E3]
set xtics 5E2


set size 1, 1.5

set lmargin 7
set bmargin 4

set multiplot

##################################################
set size 1, 0.75
set origin 0,0.75

set xlabel ""
set format x ""

set ylabel "Inverse participation ratio" offset 0.8, 0
set yrange [0 to 0.1]
set ytics 0.05


set label 1 at graph 0.1, 0.95 labstr

plot rt every 15 u 1:5 w l lw 2


##################################################
set size 1,0.9
set origin 0,0

unset label 1

set xlabel xlab
set format x "%g"

rl="Re[{/Symbol l}]"
il="Im[{/Symbol l}]"


set key top right

set ylabel rl.", ".il offset 0.8, 0
set yrange [-2 to 2]
set ytics 2
plot rt every 15 u 1:2 w l lw 2 lt 1 t rl, \
     rt every 15 u 1:3 w l lw 2 lt 3 t il
##################################################
set size 0.5, 0.6
set origin 0.45, 0.85

set grid front
unset grid
xmin=0.1
xmax=1000
set xrange [xmin to xmax]
set logscale x
set xtics 100
set mxtics 0
set xlabel "{/Symbol n}, kHz"
set format x "10^{%T}"

ymin=1E-5
ymax=1E2
set ylabel "|{/Symbol l}({/Symbol n})|^2"
set yrange [ymin to ymax]
set logscale y
set ytics 100
set mytics 0
set format y "10^{%T}"

# Bare recoil frequency
f0=2*pi*2*50
#set arrow from first f0, ymax to first f0, ymin lc rgb("red") lw 2 front


unset key
set obj rect from graph 0,0 to graph 1,1 fs solid fc rgb("white") back

plot ft u 1:2 w l lc rgb("gray") lw 4


unset multiplot

set macros
!epstopdf @fname
