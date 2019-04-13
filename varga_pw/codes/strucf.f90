subroutine structure_factor
  USE GVECTOR
  USE PSEUDOPOTENTIAL
      implicit none
      complex*16        :: stf
      integer           :: ia, is, ig,ii1,ii2,ii3
      integer, external :: iflip
      do  is=1,n_species
        do  ig=1,n_G_vector_max
          ii1=iflip(G_vector(1,ig),N_L(1))
          ii2=iflip(G_vector(2,ig),N_L(2))
          ii3=iflip(G_Vector(3,ig),N_L(3))
          stf=(0.0,0.0)
          do  ia=1,n_atom(is)
            stf = stf +ei1(ii1,ia,is)*ei2(ii2,ia,is)*ei3(ii3,ia,is)
          end do
          sfac(is,ig)=stf
        end do
      end do
end subroutine structure_factor
