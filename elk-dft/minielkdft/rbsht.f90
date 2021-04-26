subroutine rbsht(nr,nri,rfmt1,rfmt2)
  use modmain
  implicit none
  ! arguments
  integer, intent(in) :: nr,nri
  real(8), intent(in) :: rfmt1(*)
  real(8), intent(out) :: rfmt2(*)
  ! local variables
  integer nro,i
  ! transform the inner part of the muffin-tin
  call dgemm('N','N',lmmaxi,nri,lmmaxi,1.d0,rbshti,lmmaxi,rfmt1,lmmaxi,0.d0, &
   rfmt2,lmmaxi)
  ! transform the outer part of the muffin-tin
  nro=nr-nri
  i=lmmaxi*nri+1
  call dgemm('N','N',lmmaxo,nro,lmmaxo,1.d0,rbshto,lmmaxo,rfmt1(i),lmmaxo,0.d0, &
   rfmt2(i),lmmaxo)
  return
end subroutine

