#!/bin/sh
gfortran -c RoeSolver_mod.f90 
gfortran RoeSolver_mod.o Euler.f90 -o Euler -ffpe-trap=invalid,zero,overflow
#./Euler
