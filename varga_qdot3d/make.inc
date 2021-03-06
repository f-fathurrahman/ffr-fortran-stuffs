EXE_MAIN = qd3d.x

MAKE = make
AR = ar

F77 = gfortran
F77_OPTS = -O3

F90 = gfortran
F90_OPTS = -cpp -Wall -O3 # -fcheck=all # -g -fbacktrace -ffpe-trap=zero,overflow,underflow
F90_OPTS_NOWARN = -cpp -O3

CC = gcc
CC_OPTS = -std=c99

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIB_XC = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)

