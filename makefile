TARGET_EXEC := prog.exe

BUILD_DIR := ./build
SRC_DIRS := ./src
FTR_COMP := gfortran
C++_COMP := g++
FLAGS := -O3
C_RANDOM := SFMT-src-1.5.1/SFMT.c

LANGUAGES = fortran cpp

all: $(LANGUAGES)
$(LANGUAGES): %: %.exe
# run_$(LANGUAGES): %: $(LANGUAGES).out
run_cpp: cpp.out
run_fortran: fortran.out

cpp.exe: main.cpp
	$(C++_COMP) $(FLAGS) main.cpp $(C_RANDOM) -DSFMT_MEXP=1279 -o $@

fortran.exe: main.f95 subroutines.o mt19937-64.o
	$(FTR_COMP) $(FLAGS) $^ -o $@

subroutines.o: subroutines.f95 mt19937-64.o
	$(FTR_COMP) $(FLAGS) -c $< -o $@

mt19937-64.o: mt19937-64.f95
	$(FTR_COMP) $(FLAGS) -c $^ -o $@




fortran.out: fortran.exe input.dat
	./fortran.exe > $@


cpp.out: cpp.exe input.dat
	./cpp.exe > $@
	# time ./cpp.exe

plot: fortran.out cpp.out
	python3 plot_results.py

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
