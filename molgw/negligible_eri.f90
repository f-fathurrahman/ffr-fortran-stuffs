!=========================================================================
subroutine negligible_eri(tol)
  USE m_definitions
  USE m_eri
 implicit none
 real(dp),intent(in) :: tol
!=====
 integer             :: icount,ibf,jbf,kbf,lbf,jcount
 integer             :: ibuffer
 real(dp)            :: integral_ij(nbf_eri,nbf_eri)
!=====
  REAL(dp) :: eri

 icount=0
 do ibuffer=1,nsize
   if( ABS( eri_4center(ibuffer) ) < tol ) icount=icount+1
 enddo

 write(stdout,*) ' number of negligible integrals <',tol
 write(stdout,*) icount, ' / ',nsize,REAL(icount,dp)/REAL(nsize,dp)*100.0_dp,' [%]'


 do ibf=1,nbf_eri
   do jbf=1,nbf_eri
     integral_ij(ibf,jbf) = eri(ibf,jbf,ibf,jbf)
   enddo
 enddo

 write(stdout,*) 'testing Cauchy-Schwarz condition'
 icount=0
 jcount=0
 do ibf=1,nbf_eri
   do jbf=1,nbf_eri
     do kbf=1,nbf_eri
       do lbf=1,nbf_eri
         if( SQRT( integral_ij(ibf,jbf) * integral_ij(kbf,lbf) ) < tol ) icount = icount + 1
         if( ABS( eri(ibf,jbf,kbf,lbf) ) < tol ) jcount = jcount + 1
       enddo
     enddo
   enddo
 enddo
 write(stdout,*) ' number of negligible integrals <',tol
 write(stdout,*) icount, ' / ',nbf_eri**4,REAL(icount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'
 write(stdout,*) jcount, ' / ',nbf_eri**4,REAL(jcount,dp)/REAL(nbf_eri,dp)**4*100.0_dp,' [%]'


end subroutine negligible_eri
