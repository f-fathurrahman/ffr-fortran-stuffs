include make.inc

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(C_SRC:.c=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.f90
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


# Libraries
lib:
	$(MAKE) -C 01_linalg_ext
	$(MAKE) -C 02_clib
	$(MAKE) -C 10_defs
	$(MAKE) -C 10_dumpinfo
	$(MAKE) -C 11_memory_mpi
	$(MAKE) -C 11_qespresso_ext
	$(MAKE) -C 12_hide_mpi
	$(MAKE) -C 14_hidewrite
	$(MAKE) -C 16_hideleave
	$(MAKE) -C 17_libtetra_ext
	$(MAKE) -C 18_timing
	$(MAKE) -C 20_datashare
	$(MAKE) -C 21_hashfuncs
	$(MAKE) -C 27_toolbox_oop
	$(MAKE) -C 28_numeric_noabirule
	$(MAKE) -C 29_xgTools
	$(MAKE) -C 30_diago
	$(MAKE) -C 32_util
	$(MAKE) -C 41_geometry
	$(MAKE) -C 41_xc_lowlevel
	$(MAKE) -C 42_libpaw
	$(MAKE) -C 42_nlstrain
	$(MAKE) -C 42_parser
	$(MAKE) -C 43_ptgroups
	$(MAKE) -C 43_wvl_wrappers
	$(MAKE) -C 44_abitypes_defs
	$(MAKE) -C 45_geomoptim
	$(MAKE) -C 49_gw_toolbox_oop
	$(MAKE) -C 51_manage_mpi
	$(MAKE) -C 52_fft_mpi_noabirule
	$(MAKE) -C 53_ffts
	$(MAKE) -C 53_spacepar
	$(MAKE) -C 54_abiutil
	$(MAKE) -C 56_io_mpi
	$(MAKE) -C 56_mixing
	$(MAKE) -C 56_recipspace
	$(MAKE) -C 56_xc
	$(MAKE) -C 57_iopsp_parser
	$(MAKE) -C 57_iovars
	$(MAKE) -C 59_ionetcdf
	$(MAKE) -C 61_occeig
	$(MAKE) -C 62_cg_noabirule
	$(MAKE) -C 62_ctqmc
	$(MAKE) -C 62_iowfdenpot
	$(MAKE) -C 62_poisson
	$(MAKE) -C 62_wvl_wfs
	$(MAKE) -C 64_atompaw
	$(MAKE) -C 64_psp
	$(MAKE) -C 65_paw
	$(MAKE) -C 66_nonlocal
	$(MAKE) -C 66_wfs
	$(MAKE) -C 67_common
	$(MAKE) -C 68_dmft
	$(MAKE) -C 68_recursion
	$(MAKE) -C 68_rsprc
	$(MAKE) -C 69_wfdesc
	$(MAKE) -C 70_gw
	$(MAKE) -C 71_bse
	$(MAKE) -C 72_response
	$(MAKE) -C 77_ddb
	$(MAKE) -C 77_suscep
	$(MAKE) -C 78_effpot
	$(MAKE) -C 78_eph
	$(MAKE) -C 79_seqpar_mpi
	$(MAKE) -C 94_scfcv
	$(MAKE) -C 95_drive
	ar rcs libs/libABINIT.a *_*/*.o


# does not delete *.x files
clean:
	cd 01_linalg_ext; make clean
	cd 02_clib; make clean
	cd 10_defs; make clean
	cd 10_dumpinfo; make clean
	cd 11_memory_mpi; make clean
	cd 11_qespresso_ext; make clean
	cd 12_hide_mpi; make clean
	cd 14_hidewrite; make clean
	cd 16_hideleave; make clean
	cd 17_libtetra_ext; make clean
	cd 18_timing; make clean
	cd 20_datashare; make clean
	cd 21_hashfuncs; make clean
	cd 27_toolbox_oop; make clean
	cd 28_numeric_noabirule; make clean
	cd 29_xgTools; make clean
	cd 30_diago; make clean
	cd 32_util; make clean
	cd 41_geometry; make clean
	cd 41_xc_lowlevel; make clean
	cd 42_libpaw; make clean
	cd 42_nlstrain; make clean
	cd 42_parser; make clean
	cd 43_ptgroups; make clean
	cd 43_wvl_wrappers; make clean
	cd 44_abitypes_defs; make clean
	cd 45_geomoptim; make clean
	cd 49_gw_toolbox_oop; make clean
	cd 51_manage_mpi; make clean
	cd 52_fft_mpi_noabirule; make clean
	cd 53_ffts; make clean
	cd 53_spacepar; make clean
	cd 54_abiutil; make clean
	cd 56_io_mpi; make clean
	cd 56_mixing; make clean
	cd 56_recipspace; make clean
	cd 56_xc; make clean
	cd 57_iopsp_parser; make clean
	cd 57_iovars; make clean
	cd 59_ionetcdf; make clean
	cd 61_occeig; make clean
	cd 62_cg_noabirule; make clean
	cd 62_ctqmc; make clean
	cd 62_iowfdenpot; make clean
	cd 62_poisson; make clean
	cd 62_wvl_wfs; make clean
	cd 64_atompaw; make clean
	cd 64_psp; make clean
	cd 65_paw; make clean
	cd 66_nonlocal; make clean
	cd 66_wfs; make clean
	cd 67_common; make clean
	cd 68_dmft; make clean
	cd 68_recursion; make clean
	cd 68_rsprc; make clean
	cd 69_wfdesc; make clean
	cd 70_gw; make clean
	cd 71_bse; make clean
	cd 72_response; make clean
	cd 77_ddb; make clean
	cd 77_suscep; make clean
	cd 78_effpot; make clean
	cd 78_eph; make clean
	cd 79_seqpar_mpi; make clean
	cd 94_scfcv; make clean
	cd 95_drive; make clean
	rm -f libs/libABINIT.a

# also delete *.x files
cleanx:
	cd extlibs/sparskit; $(MAKE) clean
	cd extlibs/poisson_ISF; $(MAKE) clean
	cd extlibs/libxc; $(MAKE) clean
	cd extlibs/Faddeeva; $(MAKE) clean
	rm -rf *.o *.mod libmain.a libpoisson_ISF.a libsparskit.a *.x


