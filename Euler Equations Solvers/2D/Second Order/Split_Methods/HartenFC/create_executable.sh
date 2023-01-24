#!/bin/sh
gfortran -c RoeSolver_mod_2D.f90 
gfortran RoeSolver_mod_2D.o EulerHarten_2D.f90 -o EulerHarten_2D -ffpe-trap=invalid,zero,overflow
./EulerHarten_2D
