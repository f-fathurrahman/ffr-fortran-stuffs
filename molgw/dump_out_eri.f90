!=========================================================================
subroutine dump_out_eri(rcut)
 implicit none
 real(dp),intent(in) :: rcut
!=====
 character(len=50) :: filename
 integer           :: nline,iline,icurrent
 integer           :: erifile
!=====

 if(rcut < 1.0e-6_dp) then
   filename='molgw_eri.data'
 else
   filename='molgw_eri_lr.data'
 endif
 write(stdout,*) 'Dump out the ERI into file'
 write(stdout,*) 'Size of file [bytes]',REAL(nsize,dp)*dp

 if( is_iomaster ) then
   open(newunit=erifile,file=TRIM(filename),form='unformatted')
   write(erifile) nsize
   write(erifile) rcut

   nline = nsize / line_length + 1
   icurrent=0
   do iline=1,nline
     write(erifile) eri_4center(icurrent+1:MIN(nsize,icurrent+line_length+1))
     icurrent = icurrent + line_length + 1
   enddo

   close(erifile)
 endif

 write(stdout,'(a,/)') ' file written'

end subroutine dump_out_eri
