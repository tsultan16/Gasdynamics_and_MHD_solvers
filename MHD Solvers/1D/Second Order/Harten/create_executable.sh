#!/bin/sh
gfortran -c vars.f90
gfortran -c init.f90
gfortran -c RoeSolver.f90 
gfortran vars.o init.o RoeSolver.o MHD_Harten.f90 -o MHD_Harten -ffpe-trap=invalid,zero,overflow
./MHD_Harten
