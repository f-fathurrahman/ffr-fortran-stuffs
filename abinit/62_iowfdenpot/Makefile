include ../make.inc

LIBNAME = lib62_iowfdenpot.a

SRC = \
initwf.F90                    m_io_screening.F90  out1dm.F90   randac.F90  WffReadEigK.F90 \
interfaces_62_iowfdenpot.F90  m_wfk.F90           outwant.F90  rdnpw.F90   WffReadSkipK.F90


include ../make.inc.post


initwf.o : initwf.F90 interfaces_62_iowfdenpot.o m_wfk.o 
interfaces_62_iowfdenpot.o : interfaces_62_iowfdenpot.F90 
m_io_screening.o : m_io_screening.F90 
m_wfk.o : m_wfk.F90 
out1dm.o : out1dm.F90 
outwant.o : outwant.F90 
randac.o : randac.F90 
rdnpw.o : rdnpw.F90 
WffReadEigK.o : WffReadEigK.F90 
WffReadSkipK.o : WffReadSkipK.F90 
