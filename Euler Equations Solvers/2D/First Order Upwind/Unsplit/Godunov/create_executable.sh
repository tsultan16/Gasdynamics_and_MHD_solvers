#!/bin/sh
gfortran -c RoeSolver_2D.f90 
gfortran -o3 RoeSolver_2D.o EulerGodunov_2D.f90 -o EulerGodunov_2D -ffpe-trap=invalid,zero,overflow
#./EulerGodunov_2D
