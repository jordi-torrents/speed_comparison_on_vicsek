#!/bin/bash

# Compile programs
make

# Default values
v0=0.03
rho=1.0
N_reset=1000
N_steps=1000
seed=1234

input_file="times_input_file.dat"

Ls="2 4 8 16"

> times_fortran.out
> times_cpp.out
> times_c.out
> times_python.out
> Ns.out

echo "Computing different system sizes..."
for L in $Ls; do
    > $input_file
    echo $L
    echo "$L*$L*$rho" | bc >> Ns.out
    echo "$L L" >> $input_file
    echo "$v0 v0" >> $input_file
    echo "$rho rho" >> $input_file
    echo "$N_reset N_reset" >> $input_file
    echo "$N_steps N_steps" >> $input_file
    echo "$seed seed" >> $input_file
    (/usr/bin/time -f "%e"  ./fortran.exe $input_file > tmp ) 2>>\
    times_fortran.out
    (/usr/bin/time -f "%e"  ./cpp.exe $input_file > tmp ) 2>>\
    times_cpp.out
    (/usr/bin/time -f "%e"  ./c.exe $input_file > tmp ) 2>>\
    times_c.out
    (/usr/bin/time -f "%e"  python3 main.py $input_file > tmp ) 2>>\
    times_python.out
done

rm tmp
rm $input_file

echo "Plotting execution times"
python3 plot_times.py
