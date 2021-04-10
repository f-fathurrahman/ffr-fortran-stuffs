
include ../make.inc

# Elk program
SRC_elk = minielkdft.f90

SRC = modmain.f90 modmpi.f90 modomp.f90 modvars.f90 \
moddftu.f90 \
modrdm.f90 \
modphonon.f90 \
modtest.f90 \
modrandom.f90 \
modscdft.f90 \
modpw.f90 \
modtddft.f90 \
modulr.f90 \
modgw.f90 \
modjx.f90 \
modw90.f90 \
libxcifc_stub.f90 \
modxcifc.f90 \
modstore.f90 \
readinput.f90 \
mkl_init.f90 mkl_stub.f90 \
oblas_init.f90 oblas_stub.f90 \
blis_init.f90 blis_stub.f90 \
closefiles.f90 \
readspecies.f90 \
r3mv.f90 \
findprimcell.f90 \
numlist.f90 \
genspecies.f90 \
genfspecies.f90 \
r3minv.f90 \
eulerrot.f90 \
axangrot.f90 \
massnucl.f90 \
sortidx.f90 \
r3cross.f90 \
r3frac.f90 \
addlorbcnd.f90 \
atom.f90 \
splintwp.f90 \
fderiv.f90 \
spline.f90 \
splint.f90 \
wsplintp.f90 \
wsplint.f90 \
potnucl.f90 \
radnucl.f90 \
polynm.f90 \
rdirac.f90 \
rdiracint.f90 \
xc_am05.f90  xc_c_tb09.f90  xc_pbe.f90  xc_pwca.f90  xc_pzca.f90  xc_vbh.f90  xc_wc06.f90  xc_xalpha.f90 \
c_pbe.f90  c_pbe_gcor.f90  x_pbe.f90  x_wc06.f90 \
rfsht.f90 rbsht.f90 \
rfint.f90 \
rfmtlm.f90 \
rfmtftoc.f90 \
rfmtsm.f90 \
ztorfmt.f90 \
ztorflm.f90 \
rtozfmt.f90 \
rtozflm.f90 \
gradrf.f90 \
gradrfmt.f90 \
gradzfmt.f90 \
zfftifc.f90 \
clebgor.f90 \
wigner3j.f90 \
wigner3jf.f90 \
wigner6j.f90 \
factr.f90 \
factnm.f90 \
fsmooth.f90 \
gndstate.f90 \
readfermi.f90 \
init2.f90 \
timesec.f90 \
gengvc.f90 \
gencfrc.f90 \
nfftifc.f90 \
writegclq.f90 \
writeinfo.f90 writeiad.f90 \
writegeom.f90 writekpts.f90 \
writelat.f90 writestate.f90 \
writeforces.f90 writespecies.f90 \
writesym.f90 \
writedftu.f90 \
writefsm.f90 \
writemom.f90 \
writechg.f90 \
writeengy.f90 \
gengclq.f90 \
findnjcmax.f90 \
writetm2du.f90 \
genppts.f90 \
writeqpts.f90 \
vecfbz.f90 \
tm2pol.f90 \
pottm2.f90 \
dmtotm2.f90 \
rhoinit.f90 \
maginit.f90 \
rfmtctof.f90 \
rfmtpack.f90 \
rfinterp.f90 \
sbessel.f90 \
sbesseldm.f90 \
projsbf.f90 \
gentau.f90 gentaucr.f90 gentauk.f90 \
symrf.f90 \
symrvf.f90 \
symrvfmt.f90 \
symrvfir.f90 \
symrfmt.f90 \
rotrfmt.f90 \
rotrflm.f90 \
rotzflm.f90 \
rlmrot.f90 \
roteuler.f90 \
ylmrot.f90 \
ylmroty.f90 \
symrfir.f90 \
zbsht.f90 \
genwfsv.f90 \
wavefmt.f90 \
wavefcr.f90 \
getevecfv.f90 \
getevecsv.f90 \
axangsu2.f90 \
rotaxang.f90 \
findkpt.f90 \
r3mtv.f90 \
match.f90 \
genylmg.f90  genylmv.f90 \
potcoul.f90 \
potefield.f90 \
sphcrd.f90 \
zpotcoul.f90 \
genzvclmt.f90 \
zpotclmt.f90 \
potxc.f90 \
potxcir.f90 \
potxcmt.f90 \
ggair_1.f90   ggair_2b.f90    ggair_sp_2a.f90  ggamt_1.f90   ggamt_2b.f90    ggamt_sp_2a.f90 \
ggair_2a.f90  ggair_sp_1.f90  ggair_sp_2b.f90  ggamt_2a.f90  ggamt_sp_1.f90  ggamt_sp_2b.f90 \
grad2rfmt.f90 \
curlrvf.f90 \
genws.f90 \
rfirsm.f90 \
init0.f90 \
reciplat.f90 \
plotpt1d.f90 \
checkmt.f90 \
genrmesh.f90 \
gengclg.f90 \
gensfacgp.f90 \
genffacgp.f90 \
gencfun.f90 \
allatoms.f90 \
nuclei.f90 \
energynn.f90 \
sdelta.f90 \
sdelta_mp.f90 \
sdelta_sq.f90 \
sdelta_lr.f90 \
sdelta_fd.f90 \
sphcover.f90 \
genrlmv.f90 \
mtdmin.f90 \
hermite.f90 \
tm2todm.f90 \
symmetry.f90 \
gengvec.f90 \
genjlgprmt.f90 \
genshtmat.f90 \
gentimes.f90 \
readafieldt.f90 \
mixerifc.f90 mixadapt.f90 \
mixbroyden.f90 mixlinear.f90 \
sctovec.f90 \
rminv.f90 zminv.f90 \
findsym.f90 \
findsymlat.f90 \
findsymcrys.f90 \
findsymsite.f90 \
r3mm.f90 \
r3mtm.f90 \
i3mdet.f90 \
i3minv.f90 \
gridsize.f90 \
init1.f90 \
findngkmax.f90 \
gengkvec.f90 \
gauntyry.f90 gaunt.f90 \
genkpakq.f90 \
genidxlo.f90 \
potks.f90 \
poly3.f90 \
poly4i.f90 \
mixpack.f90 \
gencore.f90 \
linengy.f90 \
writelinen.f90 \
genapwlofr.f90 \
gensocfr.f90 \
genevfsv.f90 \
occupy.f90 \
putoccsv.f90 \
writeeval.f90 \
writefermi.f90 \
rhomag.f90 \
genvsig.f90 \
energy.f90 \
rfinp.f90 \
energykncr.f90 \
exxengy.f90 \
exxengyk.f90 \
zfsht.f90 \
zfmtinp.f90 \
getkmat.f90 \
zfinp.f90 \
genzrho.f90 \
rfmtinp.f90 \
rhomagsh.f90 \
rhomagk.f90 \
stheta.f90 \
charge.f90 \
findswidth.f90 \
rhocore.f90 \
moment.f90 \
gengclgq.f90 \
rhonorm.f90 \
stheta_sq.f90 \
stheta_lr.f90 \
stheta_mp.f90 \
stheta_fd.f90 \
eveqn.f90 \
putevalfv.f90 \
putevalsv.f90 \
putevecfv.f90 \
putevecsv.f90 \
chargemt.f90 \
rfpack.f90 \
findband.f90 \
eveqnit.f90 \
eveqnfv.f90 \
eveqnss.f90 \
eveqnsv.f90 \
eveqnz.f90 \
lopzflm.f90 \
eveqnfvz.f90 \
eveqnfvr.f90 \
hmlfv.f90 \
olpfv.f90 \
rschrodint.f90 \
hmlaa.f90 \
hmlalo.f90 \
hmllolo.f90 \
hmlistl.f90 \
olpaa.f90 \
olpalo.f90 \
olplolo.f90 \
olpistl.f90 \
getevalfv.f90 \
zmctmu.f90 \
olprad.f90 \
hmlrad.f90 \
genapwfr.f90 \
genlofr.f90 \
erf.f90 \
init_am_variables.f90 \
init_idx_atom_species.f90 \
init_spin_variables.f90 \
init_crystal_structure.f90 \
init_vector_field_E_A.f90 \
init_radial_meshes.f90 \
minielkdft.f90

OBJ = $(SRC:.f90=.o)
EXE = minielkdft.x

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
lib:	$(OBJ)
	ar rcs libelk.a $(OBJ)

main:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_W90) $(LIB_LPK) $(LIB_FFT) \
    $(LIB_libxc) $(LIB_SYS)

clean:
	rm -f *.o *.mod 

