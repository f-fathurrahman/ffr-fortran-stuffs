gfortran "$1" -I SRC_octopus-9.2/src/ -I SRC_octopus-9.2/src/include/ \
SRC_octopus-9.2/src/liboctopus.a \
SRC_octopus-9.2/liboct_parser/liboct_parser.a \
SRC_octopus-9.2/external_libs/qshep/libqshep.a \
SRC_octopus-9.2/external_libs/dftd3/libdftd3.a \
SRC_octopus-9.2/external_libs/bpdn/libbpdn.a \
SRC_octopus-9.2/external_libs/spglib-1.9.9/src/libsymspg.a \
-lgsl -lgslcblas -lstdc++ -lblas -llapack -lxcf90 -lxc -lfftw3
