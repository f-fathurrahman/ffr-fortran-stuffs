#!/bin/bash
FC=gfortran

OPTS="-Wall -I../"

LIBS="../libmain.a -lblas -llapack /home/efefer/mysoftwares/libint-2.4.2/lib/libint2.a -lstdc++"

basnam=`basename $1 .f90`
$FC $OPTS $1 -o $basnam.x $LIBS

