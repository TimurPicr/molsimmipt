#!/bin/bash
set -e
if [ "$#" -eq 0 ]; then echo "Usage: $0 LOG [LOG ...]" >&2; exit 1; fi

energy=""; pressure=""; temperature=""
for log in "$@"; do
    awk '$1=="Step" {e=p=t=0; for(i=1;i<=NF;i++) {
             if($i=="v_energy_peratom") e=i; if($i=="Press") p=i; if($i=="Temp") t=i }
         if(!e || !p || !t) {print "Missing energy, pressure or temperature column" > "/dev/stderr"; exit 1} data=1; next}
         data && $1~/^[0-9]+$/ {print $1*0.001, $e, $p, $t; next}
         data {exit}' "$log" > "$log.thermo"
    energy="${energy}${energy:+, }\"$log.thermo\" u 1:2 w l title \"$log\""
    pressure="${pressure}${pressure:+, }\"$log.thermo\" u 1:3 w l title \"$log\""
    temperature="${temperature}${temperature:+, }\"$log.thermo\" u 1:4 w l title \"$log\""
done

gnuplot --persist -e "set multiplot layout 3,1; set ylabel 'Atomic energy [eV]'; plot $energy; set ylabel 'Pressure [bar]'; plot $pressure; set xlabel 'Time [ps]'; set ylabel 'Temperature [K]'; plot $temperature; unset multiplot"
