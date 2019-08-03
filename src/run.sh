#/bin/bash

rm -f out
gfortran mersenne_twister.f90 precision.f03 binning.f03 MonteCarlo.f90 -o out -O1 -ggdb
#gfortran precision.f03 binning.f03 MonteCarlo.f90 mersenne_twister.f90 
./out 7172 apple .True. 10.0
