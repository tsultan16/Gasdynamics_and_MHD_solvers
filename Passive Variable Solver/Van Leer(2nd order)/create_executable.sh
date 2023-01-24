#!/bin/sh
#gfortran -c constants.f90
#gfortran -fopenmp -c RoeSolver_2D.f90
#gfortran -c velocityTracerModule.f90
#gfortran -c massTracerModule.f90
#gfortran -O2 constants.o massTracerModule.o RoeSolver_2D.o velocityTracerModule.o #EulerMH_mass.f90  -o EulerMH_mass -ffpe-trap=invalid,zero,overflow
#./EulerMH_mass

#gfortran -fopenmp constants.o RoeSolver_2D.o velocityTracerModule.o MHD_MH.f90 -o MHD_MH -#ffpe-trap=invalid,zero,overflow

#gfortran -fopenmp constants.o RoeSolver_2D.o velocityTracerModule.o massTracerModule.o MHD_MH.f90 -o MHD_MH -ffpe-trap=invalid,zero,overflow

#./MHD_MH


gfortran van_leer_pass.f90 -o pass
./pass
gnuplot "plotOutput.p"
