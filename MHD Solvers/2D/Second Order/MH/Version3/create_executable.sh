#!/bin/sh
gfortran -fopenmp -c RoeSolver_2D.f90 
gfortran -fopenmp RoeSolver_2D.o MHD_MH.f90 -o MHD_MH -ffpe-trap=invalid,zero,overflow 
#./MHD_MH
