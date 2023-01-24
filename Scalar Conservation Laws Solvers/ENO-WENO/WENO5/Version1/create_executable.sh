#!/bin/sh -l
gfortran weno5.f90 -o weno5
./weno5
gnuplot "plotOutput.p"
