#!/bin/sh
gfortran -c vars.f90
gfortran -c init.f90
gfortran -c RoeSolver.f90 
gfortran vars.o init.o RoeSolver.o MHDMH.f90 -o MHDMH -ffpe-trap=invalid,zero,overflow
./MHDMH
