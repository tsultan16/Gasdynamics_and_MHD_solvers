#!/bin/sh
gfortran -c vars.f90
gfortran -c init.f90
gfortran -c RoeSolver.f90 
gfortran vars.o init.o RoeSolver.o MHD_Balsara.f90 -o MHDROE -ffpe-trap=invalid,zero,overflow
./MHD_Balsara
