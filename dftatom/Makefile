F90 = gfortran
F90FLAGS = -Wall -O3

PROG = libmain.a

OBJS = stop_error.o \
gen_rmesh_exp.o gen_drmesh_exp.o \
integrate_trapz.o

all: $(PROG)

$(PROG): $(OBJS)
	ar rcs $@ $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

