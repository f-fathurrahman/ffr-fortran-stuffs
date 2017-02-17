#!/bin/bash
#FC=gfortran
#

FC=ifort
OPTS="-warn -no-gen-interfaces"

basnam=`basename $1 .f90`
$FC $OPTS $1 -o $basnam.x

