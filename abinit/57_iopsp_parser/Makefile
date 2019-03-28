include ../make.inc

LIBNAME = lib57_iopsp_parser.a

SRC = \
m_psxml2ab.F90 \
inpspheads.F90  \
interfaces_57_iopsp_parser.F90  pawpsxml2ab.F90  upfheader2abi.F90

include ../make.inc.post

inpspheads.o : inpspheads.F90 interfaces_57_iopsp_parser.o m_psxml2ab.o 
interfaces_57_iopsp_parser.o : interfaces_57_iopsp_parser.F90 
m_psxml2ab.o : m_psxml2ab.F90 
pawpsxml2ab.o : pawpsxml2ab.F90 
upfheader2abi.o : upfheader2abi.F90 interfaces_57_iopsp_parser.o

