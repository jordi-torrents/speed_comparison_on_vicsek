#!/bin/bash

make

# Default values
v0=0.03
rho=1.0
N_reset=10000
N_steps=10000
seed=1234

input_file="times_input_file.dat"

for L in 1 2 3 4 5
do
    echo $L > $input_file
    echo $v0 >> $input_file
    echo $rho >> $input_file
    echo $N_reset >> $input_file
    echo $N_steps >> $input_file
    echo $seed >> $input_file
    # echo ./fortran.exe $input_file
    time ./fortran.exe $input_file > tmp
done

rm tmp
rm $input_file
