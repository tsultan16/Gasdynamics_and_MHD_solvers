#!/bin/sh
gfortran -c RoeSolver_entropyfixed.f90 
gfortran RoeSolver_entropyfixed.o EulerHartenFC.f90 -o EulerHartenFC
./EulerHartenFC
