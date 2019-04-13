      subroutine form_factor

      use Gvector
      USE PSEUDOPOTENTIAL
      USE PW
      implicit none

      integer      :: is, ig, ir,ll
      real*8       :: facc, vv, rr, rr2, v, integ, g_old, g_new, arg
      real*8       :: vloc(N_PP_max),psi(N_PP_max)


! the self-energy E_self of the ionic pseudocharges

      E_self=0.0
      do is=1,n_species
        E_self=E_self+charge_pp(is)*charge_pp(is)/beta*n_atom(is)
      enddo
      E_self=E_self/sqrt(2.d0*pi)

!  form factors of the ionic pseudocharge (rhops)

      do is=1,n_species
        facc=0.25 * tpiba2*beta**2
        do ig=2,N_G_vector_max
          rhops(is,ig)=-charge_pp(is)/volume * exp(-facc*g_vector_length(ig))
        enddo
      enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! the form factors of pseudopotential (vps)
!  = local part of ps-pot plus field of gaussian pseudo-charge
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do is=1,n_species
        ll = l_loc(is)

          do ir=1,N_PP(is)
            vv=-charge_pp(is)/r(ir,is)*derf(r(ir,is)/beta)
            vloc(ir)=vion(ir,is,ll)-vv
          enddo

! g = 0 component of pseudo-pot

        do ir=1,N_PP(is)
          psi(ir)=vloc(ir)*r(ir,is)**3
        enddo
        integ=sum(psi(1:N_PP(is)))*clog(is)
        vps(is,1)=integ*4.d0*pi/volume


        g_old = -1.0


        do ig=2,N_G_vector_max
          g_new = g_vector_length(ig)
          if (g_new.ne.g_old) then
            g_old = g_new
            do ir=1,N_PP(is)
              arg=r(ir,is)*sqrt(g_new)*tpiba
              psi(ir)=vloc(ir)*sin(arg)/arg*r(ir,is)**3
            enddo
            integ=sum(psi(1:N_PP(is)))*clog(is)
          endif
          vps(is,ig)=integ*4.d0*pi/volume
        enddo
      enddo
      end


      subroutine nl_pp_form_factor
      USE GVECTOR
      USE PSEUDOPOTENTIAL
      USE PW
      implicit none
      integer is, ir, ll, i_l, m_l, ik, ig, igp, i, j
      real*8  fpibo, rr2, v, rr, dv, integ, t, arg, fac,ag
      real*8  psi(N_PP_max,N_species),psi1(N_PP_max)
      real*8  cosx,cosy,cosz,gg(3)
      integer i_m,j_m,mlmax
      integer,external  :: iflip
      integer   :: ii1,ii2,ii3
      
      do  is=1,n_species
      m_l=0
      do i_l=1,l_max(is)
      if(i_l.ne.l_loc(is)) then
              do ir=1,N_PP(is)
                rr=r(ir,is)
                rr2=rr*rr
                dv=vion(ir,is,i_l)-vion(ir,is,l_loc(is))
                psi(ir,is)=dv*p(ir,is,i_l)*rr2
                psi1(ir)=dv*p(ir,is,i_l)**2*rr 
              end do
       integ=sum(psi1(1:N_PP(is)))*clog(is)
       do j=1,2*i_l-1
          wnl(is,m_l+j)=(2.0*i_l-1.0)*4.d0*pi/volume/integ
       end do
       do ik=1,N_k_points
        do ig=1,n_g_vector(ik)

          ii1=iflip(G_vector(1,G_index(ig,ik)),N_L(1))
          ii2=iflip(G_vector(2,G_index(ig,ik)),N_L(2))
          ii3=iflip(G_vector(3,G_index(ig,ik)),N_L(3))
          gg(:) = ii1*R_lattice_vector(:,1)+ &
&                 ii2*R_lattice_vector(:,2)+ &
&                 ii3*R_lattice_vector(:,3)+k_point(:,ik)
            
          t=sqrt(gplusk(ig,ik))
          arg=t*tpiba
          if (i_l.eq.1) then
            if (t.lt.1.e-4) then
              do ir=1,N_PP(is)
                psi1(ir)=psi(ir,is)
              end do
            else
              do ir=1,N_PP(is)
                fac=arg*r(ir,is)
                fac=sin(fac)/fac
                psi1(ir)=fac*psi(ir,is)
              end do
            endif
            
            integ=sum(psi1(1:N_PP(is)))*clog(is)
            pkg_a(m_l+1,ig,is,ik)=integ
          endif
          if (i_l.eq.2) then    
            if (t.lt.1.e-4) then   
              do j=1,3
                pkg_a(m_l+j,ig,is,ik)=0.0
              end do
            else
              do ir=1,N_PP(is)
                fac=arg*r(ir,is)
                fac=(sin(fac)/fac-cos(fac))/fac
                psi1(ir)=fac*psi(ir,is)
               end do
               integ=sum(psi1(1:N_PP(is)))*clog(is)
               pkg_a(m_l+1,ig,is,ik)=integ*gg(1)/t
               pkg_a(m_l+2,ig,is,ik)=integ*gg(2)/t
               pkg_a(m_l+3,ig,is,ik)=integ*gg(3)/t
            endif
          endif
          if (i_l.eq.3) then
            if (t.lt.1.e-4) then
              do j=1,5
                pkg_a(ig,is,ik,m_l+j)=0.0
              end do
            else
            do ir=1,N_PP(is)
              ag=arg*r(ir,is)
              fac=(3.0/(ag**2.0)-1.0)*sin(ag)-3.0/ag*cos(ag)
              fac=fac/ag
              psi1(ir)=fac*psi(ir,is)
            end do
            integ=sum(psi1(1:N_PP(is)))*clog(is)
            cosx=gg(1)/t
            cosy=gg(2)/t
            cosz=gg(3)/t
            pkg_a(m_l+1,ig,is,ik)=integ*0.5*(3.0*cosz**2.0-1.0)
            pkg_a(m_l+2,ig,is,ik)=integ*sqrt(3.d0)*cosz*cosx
            pkg_a(m_l+3,ig,is,ik)=integ*sqrt(3.d0)*cosz*cosy
            pkg_a(m_l+4,ig,is,ik)=integ*sqrt(3.d0)*cosx*cosy
            pkg_a(m_l+5,ig,is,ik)=integ*sqrt(3.d0)/2.0*(cosx**2.0-cosy**2.0)                    
          endif
          endif
       end do
       end do
       m_l=m_l+2*i_l-1
       endif
       end do
       end do
 
      end
      
      subroutine nl_pp
      USE GVECTOR
      USE PSEUDOPOTENTIAL
      USE PW

      implicit none



      integer 		:: k, is, ia, mlmax, iik, ik, m, i, ii, ig, igp
      integer		:: i_lm, lm_end
      real*8     	:: e_nl
      complex*16	:: t_eigr_pkg(N_G_K_vector_max,(L_PP_max+1)**2)
      complex*16	:: sf0, tt1



      integer		:: i_c, iani, iai
      complex*16	:: ct,sf


      E_non_local=0.0

      do is=1,n_species

        do ia = 1,n_atom(is)
            do iik=1,N_k_points
              ik = iik
              lm_end=l_max(is)**2+1-2*l_loc(is)
              do i_lm=1,lm_end        
                 do ig = 1,n_g_vector(ik)
                    t_eigr_pkg(ig,i_lm)=conjg(eigr(ig,ia,is,ik))*pkg_a(i_lm,ig,is,ik)
                 enddo
                 do i=1,n_orbitals
                   ct=(0.d0,0.d0)
                   do ig=1,n_g_vector(ik)
                     ct=ct+wave_function_c(ig,i,ik)*t_eigr_pkg(ig,i_lm)
                   end do
                   fnl(i,ik,is,ia,i_lm)=ct                   
                 end do                 
                 e_nl = 0.0
                 do i=1,n_orbitals
                    sf = fnl(i,ik,is,ia,i_lm)
                    e_nl = e_nl + occupation(i,ik) *   (dble(sf)**2+aimag(sf)**2)
                 enddo
                E_non_local = E_non_local + N_sym*w_k_point(ik)*e_nl * wnl(is,i_lm)
            enddo 
            enddo 
        enddo 
      enddo 

      end

