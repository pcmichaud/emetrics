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
MAIN =  runtime/runmodel 
SRC = $(OBJ)/model.o src/runmodel.f90 
all: $(MAIN) 
	@echo 'Compiled Code !'

$(MAIN): $(SRC)
	$(FORTRAN_COMPILER) $(SRC) -J$(OBJ) -O3 -o runtime/runmodel  $(FFLAGS)  

$(OBJ)/model.o: src/model.f90
	$(FORTRAN_COMPILER)  -O3 -J$(OBJ) -c $< -o $@  $(FFLAGS) 

clean:
	rm -f $(OBJ)/*.mod $(OBJ)/*.o runtime/runmodel output/*.dat