
! Copyright (C) 2002-2009 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modmain

use m_lattice
use m_atomic
use m_mt_rad_am
use m_spin
use m_electric_vector_pot
use m_symmetry
use m_gvectors
use m_kpoints
use m_gkvectors
use m_qpoints
use m_sht
use m_density_pot_xc
use m_mixing
use m_charge_moment_current
use m_apwlo
use m_hamiltonian
use m_states
use m_core_states
use m_energy
use m_force_stress




!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
integer maxscl
! current self-consistent loop number
integer iscl
! Kohn-Sham potential convergence tolerance
real(8) epspot
! energy convergence tolerance
real(8) epsengy
! force convergence tolerance
real(8) epsforce
! stress tensor convergence tolerance
real(8) epsstress

!----------------------------------------------------------!
!     density of states, optics and response variables     !
!----------------------------------------------------------!
! number of energy intervals in the DOS/optics function plot
integer nwplot
! fine k-point grid size for integration of functions in the Brillouin zone
integer ngrkf
! smoothing level for DOS/optics function plot
integer nswplot
! energy interval for DOS/optics function plot
real(8) wplot(2)
! maximum angular momentum for the partial DOS plot
integer lmaxdos
! dosocc is .true. if the DOS is to be weighted by the occupancy
logical dosocc
! dosmsum is .true. if the partial DOS is to be summed over m
logical dosmsum
! dosssum is .true. if the partial DOS is to be summed over spin
logical dosssum
! number of optical matrix components required
integer noptcomp
! required optical matrix components
integer optcomp(3,27)
! intraband is .true. if the intraband term is to be added to the optical matrix
logical intraband
! lmirep is .true. if the (l,m) band characters should correspond to the
! irreducible representations of the site symmetries
logical lmirep
! spin-quantisation axis in Cartesian coordinates used when plotting the
! spin-resolved DOS (z-axis by default)
real(8) sqados(3)
! q-vector in lattice and Cartesian coordinates for calculating the matrix
! elements < i,k+q | exp(iq.r) | j,k >
real(8) vecql(3),vecqc(3)
! maximum initial-state energy allowed in ELNES transitions
real(8) emaxelnes
! structure factor energy window
real(8) wsfac(2)

!-------------------------------------!
!     1D/2D/3D plotting variables     !
!-------------------------------------!
! number of vertices in 1D plot
integer nvp1d
! total number of points in 1D plot
integer npp1d
! vertices in lattice coordinates for 1D plot
real(8), allocatable :: vvlp1d(:,:)
! distance to vertices in 1D plot
real(8), allocatable :: dvp1d(:)
! plot vectors in lattice coordinates for 1D plot
real(8), allocatable :: vplp1d(:,:)
! distance to points in 1D plot
real(8), allocatable :: dpp1d(:)
! corner vectors of 2D plot in lattice coordinates
real(8) vclp2d(3,0:2)
! grid sizes of 2D plot
integer np2d(2)
! corner vectors of 3D plot in lattice coordinates
real(8) vclp3d(3,0:3)
! grid sizes of 3D plot
integer np3d(3)

!----------------------------------------!
!     OEP and Hartree-Fock variables     !
!----------------------------------------!
! maximum number of core states over all species
integer ncrmax
! maximum number of OEP iterations
integer maxitoep
! OEP step size
real(8) tauoep
! magnitude of the OEP residual
real(8) resoep
! exchange potential and magnetic field
real(8), allocatable :: vxmt(:,:),vxir(:)
real(8), allocatable :: bxmt(:,:,:),bxir(:,:)
! hybrid is .true. if a hybrid functional is to be used
logical hybrid
! hybrid functional mixing coefficient
real(8) hybridc

!-------------------------------------------------------------!
!     response function and perturbation theory variables     !
!-------------------------------------------------------------!
! |G| cut-off for response functions
real(8) gmaxrf
! energy cut-off for response functions
real(8) emaxrf
! number of G-vectors for response functions
integer ngrf
! number of response function frequencies
integer nwrf
! complex response function frequencies
complex(8), allocatable :: wrf(:)
! maximum number of spherical Bessel functions on the coarse radial mesh over
! all species
integer njcmax

!-------------------------------------------------!
!     Bethe-Salpeter equation (BSE) variables     !
!-------------------------------------------------!
! number of valence and conduction states for transitions
integer nvbse,ncbse
! default number of valence and conduction states
integer nvbse0,ncbse0
! maximum number of extra valence and conduction states
integer, parameter :: maxxbse=20
! number of extra valence and conduction states
integer nvxbse,ncxbse
! extra valence and conduction states
integer istxbse(maxxbse),jstxbse(maxxbse)
! total number of transitions
integer nvcbse
! size of blocks in BSE Hamiltonian matrix
integer nbbse
! size of BSE matrix (= 2*nbbse)
integer nmbse
! index from BSE valence states to second-variational states
integer, allocatable :: istbse(:,:)
! index from BSE conduction states to second-variational states
integer, allocatable :: jstbse(:,:)
! index from BSE valence-conduction pair and k-point to location in BSE matrix
integer, allocatable :: ijkbse(:,:,:)
! BSE Hamiltonian
complex(8), allocatable :: hmlbse(:,:)
! BSE Hamiltonian eigenvalues
real(8), allocatable :: evalbse(:)
! if bsefull is .true. then the full BSE Hamiltonian is calculated, otherwise
! only the Hermitian block
logical bsefull
! if hxbse/hdbse is .true. then the exchange/direct term is included in the BSE
! Hamiltonian
logical hxbse,hdbse

!--------------------------!
!     timing variables     !
!--------------------------!
! initialisation
real(8) timeinit
! Hamiltonian and overlap matrix set up
real(8) timemat
! first-variational calculation
real(8) timefv
! second-variational calculation
real(8) timesv
! charge density calculation
real(8) timerho
! potential calculation
real(8) timepot
! force calculation
real(8) timefor

!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8), parameter :: fourpi=12.566370614359172954d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
! complex constants
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
complex(8), parameter :: zi=(0.d0,1.d0)
! array of i^l and (-i)^l values
complex(8), allocatable :: zil(:),zilc(:)
! Pauli spin matrices:
! sigma_x = ( 0  1 )   sigma_y = ( 0 -i )   sigma_z = ( 1  0 )
!           ( 1  0 )             ( i  0 )             ( 0 -1 )
complex(8) sigmat(2,2,3)
data sigmat / (0.d0,0.d0), (1.d0,0.d0), (1.d0,0.d0), (0.d0,0.d0), &
              (0.d0,0.d0), (0.d0,1.d0),(0.d0,-1.d0), (0.d0,0.d0), &
              (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0),(-1.d0,0.d0) /
! Planck constant in SI units (exact, CODATA 2018)
real(8), parameter :: h_si=6.62607015d-34
! reduced Planck constant in SI units
real(8), parameter :: hbar_si=h_si/twopi
! speed of light in SI units (exact, CODATA 2018)
real(8), parameter :: sol_si=299792458d0
! speed of light in atomic units (=1/alpha) (CODATA 2018)
real(8), parameter :: sol=137.035999084d0
! scaled speed of light
real(8) solsc
! Hartree in SI units (CODATA 2018)
real(8), parameter :: ha_si=4.3597447222071d-18
! Hartree in eV (CODATA 2018)
real(8), parameter :: ha_ev=27.211386245988d0
! Hartree in inverse meters
real(8), parameter :: ha_im=ha_si/(h_si*sol_si)
! Boltzmann constant in SI units (exact, CODATA 2018)
real(8), parameter :: kb_si=1.380649d-23
! Boltzmann constant in Hartree/kelvin
real(8), parameter :: kboltz=kb_si/ha_si
! electron charge in SI units (exact, CODATA 2018)
real(8), parameter :: e_si=1.602176634d-19
! Bohr radius in SI units (CODATA 2018)
real(8), parameter :: br_si=0.529177210903d-10
! Bohr radius in Angstroms
real(8), parameter :: br_ang=br_si*1.d10
! atomic unit of magnetic flux density in SI
real(8), parameter :: b_si=hbar_si/(e_si*br_si**2)
! atomic unit of electric field in SI
real(8), parameter :: ef_si=ha_si/(e_si*br_si)
! atomic unit of time in SI
real(8), parameter :: t_si=hbar_si/ha_si
! electron g-factor (CODATA 2018)
real(8), parameter :: gfacte=2.00231930436256d0
! electron mass in SI (CODATA 2018)
real(8), parameter :: em_si=9.1093837015d-31
! atomic mass unit in SI (CODATA 2018)
real(8), parameter :: amu_si=1.66053906660d-27
! atomic mass unit in electron masses
real(8), parameter :: amu=amu_si/em_si

!---------------------------------!
!     miscellaneous variables     !
!---------------------------------!
! code version
integer version(3)
data version / 6,3,2 /
! maximum number of tasks
integer, parameter :: maxtasks=40
! number of tasks
integer ntasks
! task array
integer tasks(maxtasks)
! current task
integer task
! tlast is .true. if the calculation is on the last self-consistent loop
logical tlast
! tstop is .true. if the STOP file exists
logical tstop
! number of self-consistent loops after which STATE.OUT is written
integer nwrite
! if wrtvars is .true. then variables are written to VARIABLES.OUT
logical wrtvars
! filename extension for files generated by gndstate
character(256) filext
! default file extension
data filext / '.OUT' /
! scratch space path
character(256) scrpath
! number of note lines
integer notelns
! notes to include in INFO.OUT
character(256), allocatable :: notes(:)

end module

