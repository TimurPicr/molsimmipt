#!/bin/bash
i=$#
args=("$@")
plots=""
let i-=1
while [ $i -ge 0 ]
do
    plots=${plots}", \"${args[$i]}.thermo\" u 3:6 w l title \"${args[$i]}\""
    awk 'BEGIN {nmax=0}
    NR <= 25 {next}
    $1 == "Step" {nmax = NR + 251} 
    NR <= nmax {print}' ${args[$i]} > ${args[$i]}.thermo
    let i-=1
done

gnuplot --persist -e "set xrange [8:31.25]; set xlabel \"Density [g/cc]\"; set ylabel \"Atomic energy [eV]\"; plot 0${plots}" 
