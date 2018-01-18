#!/bin/bash
FC=pgf90
OPTS=""

basnam=`basename $1 .cuf`
$FC $OPTS $1 -o ${basnam}_gpu.x

