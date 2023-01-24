#!/bin/sh
gfortran -c RoeSolver_2D.f90 
gfortran RoeSolver_2D.o EulerENO_2D.f90 -o EulerENO_2D -ffpe-trap=invalid,zero,overflow
#./EulerENO_2D
