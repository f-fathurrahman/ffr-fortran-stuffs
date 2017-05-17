## G95

No need for additional compiler flags

Using g95:
```
a.out --g95 images=2
```

## Intel Fortran
For use with Intel Fortran compiler, enable Intel MPI:
```
source ~/intel/impi/2017.1.132/bin64/mpivars.sh
```

Compile and link:
```
ifort -coarray hello.f90
```

No need for additional flags when running the executable.
The default number of images (processor) is the maximum number of available processors.

I only tested for shared version of Intel's coarray Fortran implementation.
The distributed version of coarray needs Intel Cluster toolkit.


## GNU Fortran

Use OpenCoarray library.



