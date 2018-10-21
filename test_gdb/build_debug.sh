#!/bin/bash
FC=gfortran
OPTS="-Wall -g"

#FC=g95
#OPTS="-Wall"

#FC=ifort
#OPTS="-warn -no-gen-interfaces"

basnam=`basename $1 .f90`
$FC $OPTS $1 -o $basnam.x

