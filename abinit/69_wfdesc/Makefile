include ../make.inc

LIBNAME = lib69_wfdesc.a

SRC = \
calc_optical_mels.F90  interfaces_69_wfdesc.F90  m_oscillators.F90  m_wfd.F90   wfd_mkrho.F90     wfd_vnlpsi.F90 \
classify_bands.F90     m_io_kss.F90              m_shirley.F90      outkss.F90  wfd_pawrhoij.F90


include ../make.inc.post


calc_optical_mels.o : calc_optical_mels.F90 m_wfd.o 
classify_bands.o : classify_bands.F90 interfaces_69_wfdesc.o m_wfd.o 
interfaces_69_wfdesc.o : interfaces_69_wfdesc.F90 m_wfd.o 
m_io_kss.o : m_io_kss.F90 
m_oscillators.o : m_oscillators.F90 
m_shirley.o : m_shirley.F90 interfaces_69_wfdesc.o m_wfd.o 
m_wfd.o : m_wfd.F90 
outkss.o : outkss.F90 m_io_kss.o 
wfd_mkrho.o : wfd_mkrho.F90 m_wfd.o 
wfd_pawrhoij.o : wfd_pawrhoij.F90 m_wfd.o 
wfd_vnlpsi.o : wfd_vnlpsi.F90 m_wfd.o 
