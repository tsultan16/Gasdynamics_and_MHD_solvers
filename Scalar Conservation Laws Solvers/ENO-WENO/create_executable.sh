#!/bin/sh -l
gfortran scalarHyperbolic.f90 -o scalarHyperbolic
./scalarHyperbolic
gnuplot "plotOutput.p"
