SUBROUTINE subdomain( x, istart, ipoints )
  IMPLICIT NONE 
  INTEGER :: istart, ipoints
  REAL(8) :: X(*)
  INTEGER :: i

  DO i = 1,ipoints
    X(istart+i) = 123.456d0
  ENDDO 

END SUBROUTINE 

SUBROUTINE sub( Npoints, X )
  USE omp_lib
  IMPLICIT NONE 
  REAL(8) :: X(*)
  INTEGER :: Npoints
  !
  INTEGER :: iam, Nt, ipoints, istart

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(X,NPOINTS)
  
  iam = omp_get_thread_num()
  Nt = omp_get_num_threads()
  
  WRITE(*,*) 'iam = ', iam

  ipoints = Npoints/Nt
  istart = iam*ipoints
  
  IF( iam == Nt-1 ) THEN 
    ipoints = Npoints - istart
  ENDIF 
  
  CALL subdomain( x, istart, ipoints )

!$OMP END PARALLEL

END SUBROUTINE 


PROGRAM main
  
  IMPLICIT NONE 
  INTEGER :: i
  INTEGER, PARAMETER :: Npoints=101
  REAL(8) :: array(Npoints)
  
  CALL sub( Npoints, array )

  WRITE(*,*) 'array = '
  DO i = 1,Npoints
    WRITE(*,*) i, array(i)
  ENDDO 

END PROGRAM 
