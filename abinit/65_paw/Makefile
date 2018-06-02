include ../make.inc

LIBNAME = lib65_paw.a

SRC = \
atomden.F90            make_fc_paw.F90       pawdensities.F90      pawpolev.F90      setnoccmmp.F90 \
calc_ubare.F90         Makefile              pawdfptenergy.F90     pawprt.F90        setrhoijpbe0.F90 \
chkpawovlp.F90         m_paw_dmft.F90        pawgrnl.F90           pawpuxinit.F90    setsymrhoij.F90 \
convert_notation.F90   m_pawfgr.F90          pawgylmg.F90          pawsushat.F90     simple_j_dia.F90 \
denfgr.F90             m_pawhr.F90           pawinit.F90           paw_symcprj.F90   smatrix_k_paw.F90 \
dsdr_k_paw.F90         m_paw_pwaves_lmn.F90  pawmkaewf.F90         pawuenergy.F90    smatrix_pawinit.F90 \
expibi.F90             m_pawpwij.F90         paw_mknewh0.F90       pawuj_det.F90     spline_paw_fncs.F90 \
fourier_interpol.F90   m_paw_slater.F90      pawmknhat.F90         pawuj_ini.F90     transgrid.F90 \
initrhoij.F90          nhatgrid.F90          pawmknhat_psipsi.F90  pawuj_red.F90     wvl_nhatgrid.F90 \
int_ang.F90            optics_paw_core.F90   pawmkrho.F90          pawuj_utils.F90 \
interfaces_65_paw.F90  optics_paw.F90        pawmkrhoij.F90        pawxenergy.F90 \
linear_optics_paw.F90  pawaccrhoij.F90       pawnabla_init.F90     qijb_kk.F90 \
make_efg_onsite.F90    pawdenpot.F90         pawnhatfr.F90         read_atomden.F90


OBJ = $(SRC:.F90=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.F90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) -c $(F77_OPTS) $<

#
# C source
#
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CC_OPTS) $<

# Targets
lib: $(OBJ)
	ar rcs $(LIBNAME) *.o; mv $(LIBNAME) ../libs; mv *.mod ../mods

# does not delete *.x files
clean:
	rm -rf *.o *.mod $(LIBNAME)



atomden.o : atomden.F90 
calc_ubare.o : calc_ubare.F90 
chkpawovlp.o : chkpawovlp.F90 
convert_notation.o : convert_notation.F90 
denfgr.o : denfgr.F90 interfaces_65_paw.o m_pawfgr.o 
dsdr_k_paw.o : dsdr_k_paw.F90 
expibi.o : expibi.F90 
fourier_interpol.o : fourier_interpol.F90 interfaces_65_paw.o m_pawfgr.o 
initrhoij.o : initrhoij.F90 
int_ang.o : int_ang.F90 
interfaces_65_paw.o : interfaces_65_paw.F90 m_paw_dmft.o m_pawfgr.o 
linear_optics_paw.o : linear_optics_paw.F90 
make_efg_onsite.o : make_efg_onsite.F90 interfaces_65_paw.o 
make_fc_paw.o : make_fc_paw.F90 
m_paw_dmft.o : m_paw_dmft.F90 
m_pawfgr.o : m_pawfgr.F90 
m_pawhr.o : m_pawhr.F90 m_paw_pwaves_lmn.o 
m_paw_pwaves_lmn.o : m_paw_pwaves_lmn.F90 
m_pawpwij.o : m_pawpwij.F90 
m_paw_slater.o : m_paw_slater.F90 
nhatgrid.o : nhatgrid.F90 
optics_paw_core.o : optics_paw_core.F90 interfaces_65_paw.o 
optics_paw.o : optics_paw.F90 interfaces_65_paw.o 
pawaccrhoij.o : pawaccrhoij.F90 
pawdenpot.o : pawdenpot.F90 interfaces_65_paw.o 
pawdensities.o : pawdensities.F90 
pawdfptenergy.o : pawdfptenergy.F90 interfaces_65_paw.o 
pawgrnl.o : pawgrnl.F90 interfaces_65_paw.o 
pawgylmg.o : pawgylmg.F90 
pawinit.o : pawinit.F90 
pawmkaewf.o : pawmkaewf.F90 interfaces_65_paw.o m_paw_pwaves_lmn.o 
paw_mknewh0.o : paw_mknewh0.F90 
pawmknhat.o : pawmknhat.F90 interfaces_65_paw.o 
pawmknhat_psipsi.o : pawmknhat_psipsi.F90 
pawmkrho.o : pawmkrho.F90 interfaces_65_paw.o m_pawfgr.o 
pawmkrhoij.o : pawmkrhoij.F90 interfaces_65_paw.o m_paw_dmft.o 
pawnabla_init.o : pawnabla_init.F90 interfaces_65_paw.o 
pawnhatfr.o : pawnhatfr.F90 
pawpolev.o : pawpolev.F90 
pawprt.o : pawprt.F90 interfaces_65_paw.o 
pawpuxinit.o : pawpuxinit.F90 interfaces_65_paw.o 
pawsushat.o : pawsushat.F90 
paw_symcprj.o : paw_symcprj.F90 
pawuenergy.o : pawuenergy.F90 
pawuj_det.o : pawuj_det.F90 interfaces_65_paw.o 
pawuj_ini.o : pawuj_ini.F90 
pawuj_red.o : pawuj_red.F90 interfaces_65_paw.o 
pawuj_utils.o : pawuj_utils.F90 interfaces_65_paw.o 
pawxenergy.o : pawxenergy.F90 
qijb_kk.o : qijb_kk.F90 
read_atomden.o : read_atomden.F90 interfaces_65_paw.o m_pawfgr.o 
setnoccmmp.o : setnoccmmp.F90 
setrhoijpbe0.o : setrhoijpbe0.F90 
setsymrhoij.o : setsymrhoij.F90 
simple_j_dia.o : simple_j_dia.F90 
smatrix_k_paw.o : smatrix_k_paw.F90 
smatrix_pawinit.o : smatrix_pawinit.F90 
spline_paw_fncs.o : spline_paw_fncs.F90 
transgrid.o : transgrid.F90 m_pawfgr.o 
wvl_nhatgrid.o : wvl_nhatgrid.F90 
