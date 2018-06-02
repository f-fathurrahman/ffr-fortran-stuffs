include ../make.inc

LIBNAME = lib41_geometry.a

SRC = \
interfaces_41_geometry.F90 \
m_crystal.F90 \
m_ab7_symmetry.F90 \
m_ewald.F90 \
m_dynmat.F90 \
bldgrp.F90             gensymspgr.F90              mkrdim.F90            symanal.F90            symmetrize_xred.F90 \
bonds_lgth_angles.F90  getptgroupma.F90            mksupercell.F90       symatm.F90             symplanes.F90 \
chiscwrt.F90           getspinrot.F90              m_supercell.F90       symaxes.F90            symptgroup.F90 \
chkdilatmx.F90         holocell.F90                polcart.F90           symbrav.F90            symredcart.F90 \
chkgrp.F90  prt_cif.F90           symcharac.F90          symrelrot.F90 \
chkorthsy.F90          ioniondist.F90              prtposcar.F90         symchk.F90             symsgcube.F90 \
chkprimit.F90          irreducible_set_pert.F90    prtspgroup.F90        symdet.F90             symsghexa.F90 \
chkrprimd.F90          littlegroup_pert.F90        ptgmadata.F90         symfind.F90            symsgmono.F90 \
dbeta.F90          randomcellpos.F90     symkpt.F90             symsgortho.F90 \
dist2.F90              remove_inversion.F90  symlatt.F90            symsgtetra.F90 \
fcart2fred.F90         shellstruct.F90       symlist_bcc.F90        symspgr.F90 \
fillcell.F90           metric.F90                  smallprim.F90         symlist_fcc.F90        xcart2xred.F90 \
fred2fcart.F90         spgdata.F90           symlist_others.F90     xred2xcart.F90 \
gensymshub4.F90        mkeuler.F90                 strconv.F90           symlist_prim.F90 \
gensymshub.F90         mkradim.F90                 stresssym.F90         symmetrize_rprimd.F90


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



