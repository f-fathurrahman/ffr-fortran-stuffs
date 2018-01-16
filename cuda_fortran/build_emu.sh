#!/bin/bash
FC=pgf90
OPTS="-Mcuda=emu"

basnam=`basename $1 .cuf`
$FC $OPTS $1 -o ${basnam}_emu.x

