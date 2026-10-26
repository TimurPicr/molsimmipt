#!/bin/bash
set -e
cd "$(dirname "$0")"
shopt -s nullglob

tmp=$(mktemp .density.XXXXXX)
trap 'rm -f "$tmp"' EXIT

for log in log.morse_pair.*; do
    [[ $log =~ ^log\.morse_pair\.([0-9]+([.][0-9]+)?)$ ]] || continue
    rcut=${BASH_REMATCH[1]}
    read -r density steps < <(awk '
        $1=="Step" {for(i=1;i<=NF;i++) if($i=="Density") col=i; data=1; next}
        data && $1~/^[0-9]+$/ {step[++n]=$1; rho[n]=$col; next}
        data {data=0}
        END {
            if(!col || !n) {print "No Density data in " FILENAME > "/dev/stderr"; exit 1}
            start=0.8*step[n]
            for(i=1;i<=n;i++) if(step[i]>=start) {sum+=rho[i]; count++}
            printf "%.8f %d\n",sum/count,step[n]
        }' "$log")
    if (( steps < 10000 )); then
        printf 'Warning: %s has only %d steps; check equilibration.\n' "$log" "$steps" >&2
    fi
    printf '%s %s\n' "$rcut" "$density" >> "$tmp"
done

[[ -s $tmp ]] || { echo 'No log.morse_pair.<cutoff> files found' >&2; exit 1; }
sort -n -k1,1 "$tmp" > density.dat
gnuplot <<'EOF'
set terminal pngcairo size 1000,650 enhanced
set output 'density.png'
set xlabel 'Morse cutoff [Angstrom]'
set ylabel 'Equilibrium density [g/cm^3]'
set grid
plot 'density.dat' using 1:2 with linespoints lw 2 pt 7 title 'Morse'
EOF
echo 'Saved density.dat and density.png'
