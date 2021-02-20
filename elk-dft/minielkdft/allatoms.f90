subroutine allatoms

use modmain, only: nspecies, xctsp, nstspmax, nstsp, nsp, nrspmax, lsp, vrsp, &
                   ksp, rhosp, occsp, rsp, evalsp, solsc, ptnucl, spzn, nrsp, spfname
use modxcifc, only: getxcdata
use modomp, only: holdthd, freethd

implicit none
logical hybrid_
integer xcspin_,xcgrad_
integer is,nthd
real(8) hybridc_
character(512) xcdescr_

! allocatable arrays
real(8), allocatable :: rwf(:,:,:)

! allocate global species charge density and potential arrays
if (allocated(rhosp)) deallocate(rhosp)
allocate(rhosp(nrspmax,nspecies))
if (allocated(vrsp)) deallocate(vrsp)
allocate(vrsp(nrspmax,nspecies))

! get the exchange-correlation functional data
call getxcdata(xctsp,xcdescr_,xcspin_,xcgrad_,hybrid_,hybridc_)
call holdthd(nspecies,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rwf) &
!$OMP NUM_THREADS(nthd)
allocate(rwf(nrspmax,2,nstspmax))
!$OMP DO
do is=1,nspecies
  write(*,*) 'Solving atom for species ', trim(spfname(is))
  call atom(solsc,ptnucl,spzn(is),nstsp(is),nsp(:,is),lsp(:,is),ksp(:,is), &
   occsp(:,is),xctsp,xcgrad_,nrsp(is),rsp(:,is),evalsp(:,is),rhosp(:,is), &
   vrsp(:,is),rwf)
end do
!$OMP END DO
deallocate(rwf)
!$OMP END PARALLEL
call freethd(nthd)

!stop 'From efefer'

return

end subroutine
!EOC

