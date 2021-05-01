#!/bin/bash

INC="-I./src/"
#LIB="./src/libelk.a ./src/fftlib.a -lblas -llapack -lfftw3"
LIB="./src/libelk.a ./src/fftlib.a -lopenblas"
#LIBXC="-L$HOME/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc"

bas=`basename $1 .f90`

# remove the previous executable
rm -vf $bas.x

mpif90 -fopenmp $INC $1 $LIB $LIBXC -o $bas.x

#gfortran -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x
echo "Test executable: $bas.x"

