#!/bin/bash

#ifort -coarray hello.f90 -L/home/efefer/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin -licaf \
#-L/home/efefer/intel/impi/2017.1.132/lib64 -lmpi_mt \
#-L/home/efefer/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin -licaf

basnam=`basename $1 .f90`

#g95 $1 -o $basnam.x

mpif90 -fcoarray=lib $1 -o $basnam.x libcaf_mpi.a

