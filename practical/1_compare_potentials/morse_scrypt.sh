#!/bin/bash
set -e
cd "$(dirname "$0")"

for rcut in 5.0 6.0 7.0 8.0 9.0 10.0; do
    lmp_serial -in lammps.input -v ptype morse -v morse_rcut "$rcut"
done
