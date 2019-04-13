
subroutine phase_factor
      use Gvector
      USE PSEUDOPOTENTIAL
      USE PW
      implicit none

      integer i, i1, i2, i3, is, ia, j, k, ik, ig, igp, iai,ii1,ii2,ii3
      real*8 sum, phi1, phi2, phi3, xmat(3,3),xmati(3,3)
      real*8 taup(3)
      integer,external :: iflip

      xmat=Lattice_vector
      call inv_r(xmat,3,xmati)
       
      do  is=1,n_species
        do  ia=1,n_atom(is)
          do i=1,3
            sum=    xmati(i,1)*position(1,ia,is)
            sum=sum+xmati(i,2)*position(2,ia,is)
            sum=sum+xmati(i,3)*position(3,ia,is)
            taup(i) = sum
          end do
          phi1 = -2.0*pi*taup(1)
          phi2 = -2.0*pi*taup(2)
          phi3 = -2.0*pi*taup(3)

! G=0
          ei1(0,ia,is)=(1.0,0.0)
          ei2(0,ia,is)=(1.0,0.0)
          ei3(0,ia,is)=(1.0,0.0)


          do  i=1,(N_L(1)+1)/2
            ei1(i,ia,is) = cmplx(dcos(i*phi1),dsin(i*phi1))
            ei1(-i,ia,is)= conjg(ei1(i,ia,is))
          end do
          do j=1,(N_L(2)+1)/2
            ei2(j,ia,is) = cmplx(dcos(j*phi2),dsin(j*phi2))
            ei2(-j,ia,is)= conjg(ei2(j,ia,is))
          end do
          
          do  k=1,(N_L(3)+1)/2
            ei3(k,ia,is) = cmplx(dcos(k*phi3),dsin(k*phi3))
            ei3(-k,ia,is) = conjg(ei3(k,ia,is))
          end do
        end do
      end do
      
      do ik=1,n_k_points
        do  is=1,n_species
          do  ia=1,n_atom(is)
            do  ig=1,n_g_vector(ik)
              igp=G_index(ig,ik)
              ii1=iflip(G_vector(1,igp),N_L(1))
              ii2=iflip(G_vector(2,igp),N_L(2))
              ii3=iflip(G_vector(3,igp),N_L(3))
              eigr(ig,ia,is,ik)=ei1(ii1,ia,is)*ei2(ii2,ia,is)*ei3(ii3,ia,is)
            end do
          end do
        end do
      end do
      
end subroutine phase_factor
