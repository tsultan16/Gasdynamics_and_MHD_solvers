#!/bin/sh -l
gfortran -c global.f90
gfortran -c weno5.f90 
gfortran global.o weno5.o driver.f90 -o driver
./driver
gnuplot "plotOutput.p"
