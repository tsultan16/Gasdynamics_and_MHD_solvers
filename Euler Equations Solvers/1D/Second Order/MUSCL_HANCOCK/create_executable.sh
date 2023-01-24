#!/bin/sh
gfortran -c RoeSolver_entropyfixed.f90 
gfortran RoeSolver_entropyfixed.o EulerMH.f90 -o EulerMH
./EulerMH
