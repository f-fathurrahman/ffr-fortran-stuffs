include ../make.inc

LIBNAME = lib71_bse.a

SRC = \
exc_build_block.F90  exc_plot.F90           m_bse_io.F90       m_exc_spectra.F90  m_hexc.F90  sumrule.F90 \
exc_build_ham.F90    interfaces_71_bse.F90  m_exc_diago.F90    m_haydock.F90      setup_bse.F90 \
exc_den.F90          m_bseinterp.F90        m_exc_itdiago.F90  m_haydock_io.F90   setup_bse_interp.F90


include ../make.inc.post


exc_build_block.o : exc_build_block.F90 m_bse_io.o 
exc_build_ham.o : exc_build_ham.F90 interfaces_71_bse.o 
exc_den.o : exc_den.F90 
exc_plot.o : exc_plot.F90 m_bse_io.o 
interfaces_71_bse.o : interfaces_71_bse.F90 
m_bseinterp.o : m_bseinterp.F90 m_haydock_io.o 
m_bse_io.o : m_bse_io.F90 
m_exc_diago.o : m_exc_diago.F90 m_exc_spectra.o m_bse_io.o m_exc_itdiago.o 
m_exc_itdiago.o : m_exc_itdiago.F90 m_bse_io.o 
m_exc_spectra.o : m_exc_spectra.F90 interfaces_71_bse.o 
m_haydock.o : m_haydock.F90 m_exc_spectra.o m_hexc.o m_bse_io.o m_haydock_io.o 
m_haydock_io.o : m_haydock_io.F90 
m_hexc.o : m_hexc.F90 m_bseinterp.o m_bse_io.o m_haydock_io.o 
setup_bse.o : setup_bse.F90 
setup_bse_interp.o : setup_bse_interp.F90 
sumrule.o : sumrule.F90 
