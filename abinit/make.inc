MAKE = make
AR = ar

F77 = gfortran
F77_OPTS = -O0 -g

F90 = gfortran
F90_OPTS = -cpp -Wall -O0 -g -I../mods -I../incs -I./ -DHAVE_CONFIG_H -ffree-line-length-none -ffree-form
F90_OPTS_NOWARN = -cpp -O0 -g

CC = gcc
CC_OPTS = -DHAVE_CONFIG_H -I. -I../incs -O0 -g

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIB_XC = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)

