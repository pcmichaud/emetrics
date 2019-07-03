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
OBJ = obj
MAIN =  runchain 
SRC = sim.f90 
all: $(MAIN) 
	@echo 'Compiled Code !'

$(MAIN): $(SRC)
	$(FORTRAN_COMPILER) $(SRC) dcdflib.a  -O3 -o runchain  $(FFLAGS)  
clean:
	rm -f runchain 
