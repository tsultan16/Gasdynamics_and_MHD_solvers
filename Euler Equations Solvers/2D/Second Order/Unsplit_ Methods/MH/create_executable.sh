#!/bin/sh
gfortran -c RoeSolver_2D.f90 
gfortran -o3 RoeSolver_2D.o EulerMH_2D.f90 -o EulerMH_2D -ffpe-trap=invalid,zero,overflow
#./EulerMH_2D
