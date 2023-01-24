#!/bin/sh
gfortran -c RoeSolver_2D.f90
gfortran -c velocityTracerModule.f90
gfortran -O2 velocityTracerModule.o RoeSolver_2D.o EulerMH_vel.f90  -o EulerMH_vel
#./EulerMH_vel

