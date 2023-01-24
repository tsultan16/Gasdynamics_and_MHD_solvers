#!/bin/sh
gfortran -c global.f90 
gfortran -c init.f90
gfortran -c wenoflux.f90
gfortran global.o init.o wenoflux.o  WENO5.f90 -o WENO5 -ffpe-trap=invalid,zero,overflow 
./WENO5
./makePlots.sh
