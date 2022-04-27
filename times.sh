#!/bin/bash

# Compile programs
make

# Default values
v0=0.03
rho=1.0
N_reset=1000
N_steps=1000
seed=1234

input_file="input_file.tmp"

Ls="2 3 4 5 6 7 8 9 10"

echo "fortran" > times_fortran.tmp
echo "cpp" > times_cpp.tmp
echo "c" > times_c.tmp
echo "python" > times_python.tmp
echo "N" > Ns.tmp
echo "L" > Ls.tmp

echo "Computing different system sizes..."
for L in $Ls; do
    > $input_file
    echo $L
    echo $L >> Ls.tmp
    echo "$L*$L*$rho" | bc >> Ns.tmp
    echo "$L L" >> $input_file
    echo "$v0 v0" >> $input_file
    echo "$rho rho" >> $input_file
    echo "$N_reset N_reset" >> $input_file
    echo "$N_steps N_steps" >> $input_file
    echo "$seed seed" >> $input_file
    (/usr/bin/time -f "%e"  ./fortran.exe $input_file > /dev/null ) 2>>\
    times_fortran.tmp
    (/usr/bin/time -f "%e"  ./cpp.exe $input_file > /dev/null ) 2>>\
    times_cpp.tmp
    (/usr/bin/time -f "%e"  ./c.exe $input_file > /dev/null ) 2>>\
    times_c.tmp
    (/usr/bin/time -f "%e"  python3 main.py $input_file > /dev/null ) 2>>\
    times_python.tmp
done

paste Ls.tmp Ns.tmp times_fortran.tmp times_cpp.tmp times_c.tmp times_python.tmp |
column -t > times.dat

# paste -d "," Ls.tmp Ns.tmp times_fortran.tmp times_cpp.tmp times_c.tmp times_python.tmp |
# sed 's/,/:,/g' |
# column -t -s: |
# sed 's/ ,/,/g' > times.dat

rm *tmp

echo "Plotting execution times"
python3 plot_times.py times.dat
