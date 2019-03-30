#!/bin/bash
FC=gfortran
OPTS="-Wall -O0 -g -cpp -DHAVE_CONFIG_H"

INC="-I../mods -I../incs"
LIBS="../libs/libABINIT.a /home/efefer/mysoftwares/libxc-3.0.0/lib/libxcf90.a
/home/efefer/mysoftwares/libxc-3.0.0/lib/libxc.a -lblas -llapack"

basnam=`basename $1 .f90`
$FC $OPTS $INC $1 -o $basnam.x $LIBS

