#!/bin/bash

INC="-I./minielkdft"
LIB="./minielkdft/libelk.a ./minielkdft/fftlib.a -lopenblas"

bas=`basename $1 .f90`

# remove the previous executable
rm -vf $bas.x

gfortran -fopenmp $INC $1 $LIB $LIBXC -o $bas.x

#gfortran -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x
echo "Test executable: $bas.x"

