#!/bin/sh -l
gfortran weno5adv.f90 -o weno5adv
gfortran vanleeradv.f90 -o vanleeradv
./weno5adv
./vanleeradv
gnuplot "plotErr.p"
