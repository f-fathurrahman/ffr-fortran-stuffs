! Calculates the radial overlap integrals of the APW and local-orbital basis
! functions. In other words, for atom $\alpha$, it computes integrals of the
! form
! $$ o^{\alpha}_{qp}=\int_0^{R_i}u^{\alpha}_{q;l_p}(r)v^{\alpha}_p(r)r^2dr $$
! and
! $$ o^{\alpha}_{pp'}=\int_0^{R_i}v^{\alpha}_p(r)v^{\alpha}_{p'}(r)r^2dr,
! \quad l_p=l_{p'} $$
! where $u^{\alpha}_{q;l}$ is the $q$th APW radial function for angular
! momentum $l$; and $v^{\alpha}_p$ is the $p$th local-orbital radial function
! and has angular momentum $l_p$.
SUBROUTINE olprad()
  USE m_atomic, ONLY: natmtot, idxis
  USE m_mt_rad_am, ONLY: nrmtmax, nrmt, wrmt
  USE m_apwlo, ONLY: lofr, nlorb, lorbl, apword, apwfr
  USE m_hamiltonian, ONLY: ololo, oalo
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ias,nr
  INTEGER :: ilo,jlo,l,io
  ! automatic arrays
  REAL(8) :: fr(nrmtmax)

  DO ias = 1,natmtot
    is = idxis(ias)
    nr = nrmt(is)
    !-------------------------------------!
    !     APW-local-orbital integrals     !
    !-------------------------------------!
    DO ilo = 1,nlorb(is)
      l = lorbl(ilo,is)
      DO io = 1,apword(l,is)
        fr(1:nr) = apwfr(1:nr,1,io,l,ias)*lofr(1:nr,1,ilo,ias)
        oalo(io,ilo,ias) = dot_product(wrmt(1:nr,is),fr(1:nr))
      ENDDO 
    ENDDO 
    !-----------------------------------------------!
    !     local-orbital-local-orbital integrals     !
    !-----------------------------------------------!
    DO ilo = 1,nlorb(is)
      l = lorbl(ilo,is)
      DO jlo = 1,nlorb(is)
        IF( lorbl(jlo,is) == l ) THEN 
          fr(1:nr) = lofr(1:nr,1,ilo,ias)*lofr(1:nr,1,jlo,ias)
          ololo(ilo,jlo,ias) = dot_product(wrmt(1:nr,is),fr(1:nr))
        ENDIF 
      ENDDO 
    ENDDO 
  ENDDO 

  RETURN 
END SUBROUTINE 


