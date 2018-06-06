MAKE = make
AR = ar

F77 = gfortran
F77_OPTS = -O3

F90 = gfortran
F90_OPTS = -cpp -Wall -O3 -ffree-line-length-none
F90_OPTS_NOWARN = -cpp -O3

LIBINT = /home/efefer/mysoftwares/libint-2.4.2

CXX = clang++
CXX_OPTS = -std=c++11 -O2 -I$(LIBINT)/include/ -I$(LIBINT)/include/libint2

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIB_XC = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)

