#!/bin/sh
gfortran -c RoeSolver_2D.f90 
gfortran RoeSolver_2D.o MHD_Roe.f90 -o MHD_Roe -ffpe-trap=invalid,zero,overflow
./MHD_Roe
