TARGET_EXEC := prog.exe

BUILD_DIR := ./build
SRC_DIRS := ./src
FTR_COMP := gfortran
C++_COMP := g++
PYTHON_PATH := /usr/bin/python3
C_COMP := gcc
FLAGS := -O3 -msse2 -DHAVE_SSE2
C_RANDOM := SFMT-src-1.5.1/SFMT.c
input_file="input_file.tmp"
# DEFAULT INPUT VALUES
v0=0.03
rho=1.0
N_reset=1000
N_steps=1000
seed=1234



LANGUAGES = fortran cpp c python
all: fortran.exe cpp.exe c.exe numba_funcs.so
# all: $(LANGUAGES)
# $(LANGUAGES): %: %.exe
# run_$(LANGUAGES): %: $(LANGUAGES).out
run_cpp: cpp.out
run_c: c.out
run_fortran: fortran.out
run_python: python.out


numba_funcs.so: numba_functions.py
	rm -f numba_funcs*
	$(PYTHON_PATH) $^
	mv numba_funcs* $@

c.exe: main.c
	$(C_COMP) $(FLAGS) main.c $(C_RANDOM) -DSFMT_MEXP=19937 -o $@ -lm

cpp.exe: main.cpp
	$(C++_COMP) $(FLAGS) main.cpp $(C_RANDOM) -DSFMT_MEXP=19937 -o $@

fortran.exe: main.f95 subroutines.o mt19937-64.o
	$(FTR_COMP) $(FLAGS) $^ -o $@

subroutines.o: subroutines.f95 mt19937-64.o
	$(FTR_COMP) $(FLAGS) -c $< -o $@

mt19937-64.o: mt19937-64.f95
	$(FTR_COMP) $(FLAGS) -c $^ -o $@




fortran.out: fortran.exe input.dat
	/usr/bin/time -f "%e" ./$^ > $@

c.out: c.exe input.dat
	/usr/bin/time -f "%e" ./$^ > $@

cpp.out: cpp.exe input.dat
	/usr/bin/time -f "%e" ./$^ > $@


python.out: main.py numba_funcs.so input.dat
	/usr/bin/time -f "%e" $(PYTHON_PATH) main.py input.dat > $@

plot: fortran.out cpp.out c.out python.out
	$(PYTHON_PATH) plot_results.py

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
	rm -f *.so
	rm -f *.times

fortran.times: fortran.exe Ls.dat
	echo "fortran" > tmp

	while read L; do\
		> $(input_file);\
		echo $$L;\
		echo "$$L L" >> $(input_file);\
		echo "$(v0) v0" >> $(input_file);\
		echo "$(rho) rho" >> $(input_file);\
		echo "$(N_reset) N_reset" >> $(input_file);\
		echo "$(N_steps) N_steps" >> $(input_file);\
		echo "$(seed) seed" >> $(input_file);\
		(/usr/bin/time -f "%e"  ./fortran.exe $(input_file) > /dev/null ) 2>> tmp;\
	done <Ls.dat
	mv tmp $@
	rm $(input_file)

cpp.times: cpp.exe Ls.dat
	echo "cpp" > tmp

	while read L; do\
		> $(input_file);\
		echo $$L;\
		echo "$$L L" >> $(input_file);\
		echo "$(v0) v0" >> $(input_file);\
		echo "$(rho) rho" >> $(input_file);\
		echo "$(N_reset) N_reset" >> $(input_file);\
		echo "$(N_steps) N_steps" >> $(input_file);\
		echo "$(seed) seed" >> $(input_file);\
		(/usr/bin/time -f "%e"  ./cpp.exe $(input_file) > /dev/null ) 2>> tmp;\
	done <Ls.dat
	mv tmp $@
	rm $(input_file)

c.times: c.exe Ls.dat
	echo "c" > tmp
	Ls="$(<Ls.dat)"
	while read L; do\
		> $(input_file);\
		echo $$L;\
		echo "$$L L" >> $(input_file);\
		echo "$(v0) v0" >> $(input_file);\
		echo "$(rho) rho" >> $(input_file);\
		echo "$(N_reset) N_reset" >> $(input_file);\
		echo "$(N_steps) N_steps" >> $(input_file);\
		echo "$(seed) seed" >> $(input_file);\
		(/usr/bin/time -f "%e"  ./c.exe $(input_file) > /dev/null ) 2>> tmp;\
	done <Ls.dat
	mv tmp $@
	rm $(input_file)

python.times: numba_funcs.so main.py Ls.dat
	echo "python" > tmp
	Ls="$(<Ls.dat)"
	while read L; do\
		> $(input_file);\
		echo $$L;\
		echo "$$L L" >> $(input_file);\
		echo "$(v0) v0" >> $(input_file);\
		echo "$(rho) rho" >> $(input_file);\
		echo "$(N_reset) N_reset" >> $(input_file);\
		echo "$(N_steps) N_steps" >> $(input_file);\
		echo "$(seed) seed" >> $(input_file);\
		(/usr/bin/time -f "%e"  $(PYTHON_PATH) main.py $(input_file) > /dev/null ) 2>> tmp;\
	done <Ls.dat
	mv tmp $@
	rm $(input_file)

Ns.times: Ls.dat
	echo "L N" > $@
	Ls="$(<Ls.dat)"
	while read L; do\
		echo -n "$$L " >> $@;\
		echo "$$L*$$L*$(rho)" | bc >> $@;\
	done <Ls.dat

times.dat:Ns.times fortran.times c.times cpp.times python.times
	paste $^ | column -t > times.dat

times: times.dat
	echo "Plotting execution times"
	$(PYTHON_PATH) plot_times.py $^

time: times
