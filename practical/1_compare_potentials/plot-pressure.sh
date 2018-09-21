#!/bin/bash
i=$#
args=("$@")
plots=""
let i-=1
while [ $i -ge 0 ]
do
    plots=${plots}"\"${args[$i]}.thermo\" u 3:(\$7 / 1e4) w l title \"${args[$i]}\", "
    awk 'BEGIN {nmax=0}
    NR <= 45 {next}
    $1 == "Step" {nmax = NR + 251} 
    NR <= nmax {print}' ${args[$i]} > ${args[$i]}.thermo
    let i-=1
done

gnuplot --persist -e "set termopt enhanced; set xrange [9:31.25]; set yrange [-25:100]; set xlabel \"Atomic volume [A^3]\"; set ylabel \"Pressure [GPa]\"; plot ${plots}0"

let i=$#-1
while [ $i -ge 0 ]
do
    rm ${args[$i]}.thermo
    let i-=1
done


