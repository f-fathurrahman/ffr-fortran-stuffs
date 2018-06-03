#!/bin/bash
FC=gfortran
OPTS="-Wall -cpp -DHAVE_CONFIG_H"

INC="-I../mods -I../incs"
LIBS="../libs/libABINIT.a -lblas -llapack"

basnam=`basename $1 .f90`
$FC $OPTS $INC $1 -o $basnam.x $LIBS

