# ============================================================================
# Name        : Makefile
# Author      : Pierre-Carl Michaud
# Version     : 0.0
# Copyright   : Your copyright notice
# Description : Simple example
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran
FFLAGS = -fbounds-check

MAIN =  runmodel 
SRC = model.o runmodel.f90 
all: $(MAIN) 
	@echo 'Compiled Code !'

$(MAIN): $(SRC)
	$(FORTRAN_COMPILER) $(SRC) -O3 -o runmodel  $(FFLAGS)  

model.o: model.f90
	$(FORTRAN_COMPILER)  -O3 -c $< -o $@  $(FFLAGS) 

clean:
	rm -f *.mod runmodel *.o *.dat