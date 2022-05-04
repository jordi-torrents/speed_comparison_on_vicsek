FTR_COMP := gfortran
C++_COMP := g++
PYTHON_PATH := python3
C_COMP := gcc
FLAGS := -O3 -msse2 -DHAVE_SSE2
C_RAND := utils/SFMT-src-1.5.1/SFMT.c
FORTRAN_RAND := utils/mt19937-64.f95
SRC := ./src
BIN := ./bin
OBJ := ./obj
TMP := ./tmp
dir_guard=@mkdir -p $(@D)

input_file="input_file.tmp"
# DEFAULT INPUT VALUES
v0=0.03
rho=1.0
N_reset=1000
N_steps=1000
seed=1234




LANGUAGES = fortran cpp c python
all: $(BIN)/c.exe $(BIN)/fortran.exe $(BIN)/cpp.exe $(BIN)/numba_funcs.so
# all: $(LANGUAGES)
# $(LANGUAGES): %: %.exe
# run_$(LANGUAGES): %: $(LANGUAGES).out
run_cpp: cpp.out
run_c: c.out
run_fortran: fortran.out
run_python: python.out


$(BIN)/numba_funcs.so: $(SRC)/numba_functions.py
	$(dir_guard)
	rm -f numba_funcs*
	$(PYTHON_PATH) $^
	mv $(SRC)/numba_funcs* $@

$(BIN)/c.exe: $(SRC)/main.c
	$(dir_guard)
	$(C_COMP) $(FLAGS) $^ $(C_RAND) -DSFMT_MEXP=19937 -o $@ -lm

$(BIN)/cpp.exe: $(SRC)/main.cpp
	$(dir_guard)
	$(C++_COMP) $(FLAGS) $^ $(C_RAND) -DSFMT_MEXP=19937 -o $@

$(BIN)/fortran.exe: $(SRC)/main.f95 $(OBJ)/subroutines.o $(OBJ)/mt19937-64.o
	$(dir_guard)
	$(FTR_COMP) $(FLAGS) -J$(OBJ) -o $@ $^
$(OBJ)/subroutines.o: $(SRC)/subroutines.f95 $(OBJ)/mt19937-64.o
	$(FTR_COMP) $(FLAGS) -J$(OBJ) -c $(SRC)/subroutines.f95 $(OBJ)/mt19937-64.o -o $@
$(OBJ)/mt19937-64.o: $(FORTRAN_RAND)
	$(dir_guard)
	$(FTR_COMP) $(FLAGS) -J$(OBJ) -c $(FORTRAN_RAND) -o $@


$(TMP)/Fortran.out: $(BIN)/fortran.exe input.dat
	$(dir_guard)
	$^ > $@

$(TMP)/C.out: $(BIN)/c.exe input.dat
	$(dir_guard)
	$^ > $@

$(TMP)/C++.out: $(BIN)/cpp.exe input.dat
	$(dir_guard)
	$^ > $@

$(TMP)/Python.out: $(SRC)/main.py $(SRC)/python_functions.py input.dat
	$(dir_guard)
	$(PYTHON_PATH) $< input.dat > $@

$(TMP)/Numba.out: $(SRC)/main.py $(BIN)/numba_funcs.so input.dat
	$(dir_guard)
	$(PYTHON_PATH) $< input.dat > $@

compare_simulation_results: results.png

results.png: $(TMP)/Fortran.out $(TMP)/C.out $(TMP)/C++.out $(TMP)/Python.out $(TMP)/Numba.out
	$(PYTHON_PATH) plot_results.py $^

clean:
	rm -f *.o *.mod *.exe *.so *.times *.tmp
	rm -fr bin



$(TMP)/Fortran.times: $(BIN)/fortran.exe Ls.dat
	@echo "Fortran" > $@
	@while read L; do\
		echo -n "\r\tComputing Fortran with size $$L";\
		echo "$$L L\n$(v0) v0\n$(rho) rho\n$(N_reset) N_reset\n$(N_steps) N_steps\n$(seed) seed" > $(input_file);\
		(/usr/bin/time -f "%e" $< $(input_file) > /dev/null ) 2>> $@;\
	done <Ls.dat
	@echo

$(TMP)/C++.times: $(BIN)/cpp.exe Ls.dat
	@echo "C++" > $@
	@while read L; do\
		echo -n "\r\tComputing C++ with size $$L";\
		echo "$$L L\n$(v0) v0\n$(rho) rho\n$(N_reset) N_reset\n$(N_steps) N_steps\n$(seed) seed" > $(input_file);\
		(/usr/bin/time -f "%e" $< $(input_file) > /dev/null ) 2>> $@;\
	done <Ls.dat
	@echo

$(TMP)/C.times: $(BIN)/c.exe Ls.dat
	@echo "C" > $@
	@while read L; do\
		echo -n "\r\tComputing C with size $$L";\
		echo "$$L L\n$(v0) v0\n$(rho) rho\n$(N_reset) N_reset\n$(N_steps) N_steps\n$(seed) seed" > $(input_file);\
		(/usr/bin/time -f "%e" $< $(input_file) > /dev/null ) 2>> $@;\
	done <Ls.dat
	@echo

$(TMP)/Python.times: $(SRC)/main.py Ls.dat
	@echo "Python" > $@
	@while read L; do\
		echo -n "\r\tComputing Python with size $$L";\
		echo "$$L L\n$(v0) v0\n$(rho) rho\n$(N_reset) N_reset\n$(N_steps) N_steps\n$(seed) seed" > $(input_file);\
		(/usr/bin/time -f "%e"  $(PYTHON_PATH) $< $(input_file) > /dev/null ) 2>> $@;\
	done <Ls.dat
	@echo

$(TMP)/Numba.times: $(BIN)/numba_funcs.so $(SRC)/main.py Ls.dat
	@echo "Numba" > $@
	@while read L; do\
		echo -n "\r\tComputing Python with size $$L";\
		echo "$$L L\n$(v0) v0\n$(rho) rho\n$(N_reset) N_reset\n$(N_steps) N_steps\n$(seed) seed" > $(input_file);\
		(/usr/bin/time -f "%e"  $(PYTHON_PATH) main.py $(input_file) numba > /dev/null ) 2>> $@;\
	done <Ls.dat
	@echo

$(TMP)/Ns.times: Ls.dat
	@echo "L N" > $@
	@while read L; do\
		echo -n "$$L " >> $@;\
		echo "$$L*$$L*$(rho)" | bc >> $@;\
	done <Ls.dat

times.dat: $(TMP)/Ns.times $(TMP)/Fortran.times $(TMP)/C.times $(TMP)/C++.times $(TMP)/Python.times $(TMP)/Numba.times
	paste $^ | column -t > times.dat

times.png: times.dat
	echo "Plotting execution times"
	$(PYTHON_PATH) plot_times.py $^

time: times.png

times: times.png




# for i in $(seq 0 1023 | tqdm --total 1024); do
#   do_something
# done

# or

# for i in $(seq 0 1023); do
#   do_something
#   echo $i
# done | tqdm --total 1024 >> /dev/null
