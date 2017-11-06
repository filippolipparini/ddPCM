!-------------------------------------------------------------------------------------
! Recall :
!
!   A_eps Phi_eps = A_oo Phi     <==  PCM
!        
!         L sigma = Phi_eps      <==  COSMO
!      
! and :
!
!   E_s = 1/2 f(eps) < Psi , sigma >  ;  F_i = - dE_s / dr_i  ;  (A_eps L)^T s = Psi
!
! Thus :
!
!   E_s' = ... < Psi           , sigma' >
!
!        = ... < (A_eps L)^* s , sigma' >
!
!        = ... <             s , A_eps L sigma' >
!
! Since A_eps' is independent of eps, let's just set A' = A_oo' = A_eps' , so that :
!
!   F = ... ( - < s , A' ( Phi - Phi_eps ) > - < A_oo^T s , Phi' > + < A_eps^T s , L' sigma > )
!
!-------------------------------------------------------------------------------------
!
subroutine compute_forces( Phi, charge, Psi, sigma, Phi_eps, f )
!
      use ddcosmo , only : zero, ngrid, nsph, nbasis, zero, lmax, intrhs,     &
                           basis, fdoka, fdokb, iprint, ncav, ui, ccav, csph, &
                           w, fdoga, eps, sprod, prtsph, wghpot
!      
      implicit none
      real*8, dimension(ncav),        intent(in)  :: Phi
      real*8, dimension(       nsph), intent(in)  :: charge
      real*8, dimension(nbasis,nsph), intent(in)  :: Psi
      real*8, dimension(nbasis,nsph), intent(in)  :: sigma
      real*8, dimension(nbasis,nsph), intent(in)  :: Phi_eps
      real*8, dimension(3,nsph),      intent(out) :: f
!
      real*8, dimension(ngrid,nsph) :: xi,g
      real*8, dimension(nbasis,nsph) :: w_lm, s, y, z
      real*8, dimension(ngrid) :: x
      real*8, dimension(nbasis) :: xlm, basloc, vplm
      real*8, dimension(3,nbasis) :: dbsloc
      real*8, dimension(lmax+1) :: vcos, vsin
      real*8 :: rvoid,e0,xx(1)
      real*8 :: ef(3,ncav),zeta(ncav),phiexp(ngrid,nsph)
!
      integer :: isph, jsph, icomp, n, i, c1, c2, cr, ii
      logical, parameter :: star=.true.
      real*8, parameter :: tokcal=627.509469d0
!
!-------------------------------------------------------------------------------------
!

!!!      call forces_pcm( Phi, Charge, Psi, Sigma, Phi_eps, F )
!!!      return

!     initialize
      f(:,:) = zero
!
!
!     STEP 1 : solve adjoint problem (A_eps L)^T s = Psi
!     --------------------------------------------------
!   
!     solve L^T y = Psi     
      call cosmo( .true., .false., xx, xx, psi, y, rvoid )
!
!     solve A_eps^T s = y
      call pcm( .true., .false., .true., xx, y, s )
!!!      call prtsph('adjoint pcm solution - new:', nsph, 0, s)
!!!      call ADJpcm( y, s )
!!!      call prtsph('adjoint pcm solution - old:', nsph, 0, s)
!
!
!     STEP 2 : compute f = - < s , A' ( Phi - Phi_eps ) >
!     ---------------------------------------------------
!
!     initialize the timer
      call system_clock( count_rate = cr )
      call system_clock( count = c1 )
!
!     weight potential
      call wghpot( Phi, g )
!
!     compute SH expansion of Phi on i-sphere
      do isph = 1,nsph
!      
        call intrhs( isph, g(:,isph), w_lm(:,isph) )
!        
      enddo

!!!      g = zero ; ii = 0
!!!!      
!!!      do isph = 1,nsph
!!!        do n = 1,ngrid
!!!!        
!!!          if ( ui(n,isph).gt.zero ) then
!!!!
!!!            ii = ii + 1
!!!!
!!!            g(n,isph) = phi(ii)
!!!!            
!!!          endif
!!!!          
!!!        enddo
!!!      enddo
!!!      do isph = 1,nsph
!!!!      
!!!        call intrhs( isph, g(:,isph), w_lm(:,isph) )
!!!!        
!!!      enddo


!
!     compute w = Phi - Phi_eps
      w_lm(:,:) = w_lm(:,:) - Phi_eps(:,:)
!
!     contract
      do isph = 1,nsph 
!      
!!!        call service_routine1_new( s, w_lm , isph, f(1:3,isph) )
        call contract_dA( s, w_lm , isph, f(1:3,isph) )
!
      enddo
!
!     flip sign 
      f = -f
!
!
!     STEP 4 : compute f = f + < y , L' sigma >
!     -----------------------------------------
!
!     expand y
      xi(:,:) = zero
!      
      do isph = 1,nsph
        do n = 1,ngrid
!        
!         compute xi
          xi(n,isph) = dot_product( y(:,isph), basis(:,n) )
!          
        enddo
      enddo
!
!     contract
      do isph = 1, nsph
!
!       accumulate f += K_a contribution to < y , L' sigma >
        call fdoka( isph, sigma, xi(:,isph), basloc, dbsloc, vplm, vcos, vsin, f(:,isph) ) 
!       accumulate f += K_b contribution to < y , L' sigma >
        call fdokb( isph, sigma, xi,         basloc, dbsloc, vplm, vcos, vsin, f(:,isph) ) 
!
      enddo
!
!
!     STEP 3 : compute f = f - < A_oo^T s , Phi' >
!     --------------------------------------------
!
!     initialize
      z(:,:) = zero ; xlm(:) = zero ; x(:) = zero ; basloc(:) = zero
      vplm(:) = zero ; vcos(:) = zero ; vsin(:) = zero
!      
!     compute z = A_oo^T s
      do isph = 1,nsph
!      
        call ADJvec( isph, zero, s(:,:), z(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!        
      enddo
!      
!     expand z
      xi(:,:) = zero
!      
      do isph = 1,nsph
        do n = 1,ngrid
!        
!         compute xi
          xi(n,isph) = dot_product( z(:,isph), basis(:,n) )
!          
        enddo
      enddo
!      
!     initialize
      ii = 0 ; phiexp = zero
!
!     loop over atoms
      do isph = 1, nsph
!
!       loop over integration points
        do n = 1, ngrid
!
!         positive contribution from integration point
          if ( ui(n,isph).gt.zero ) then
!                  
!           advance index
            ii = ii + 1
!
!           expand potential Phi_n^j on a sphere-by-sphere basis [ needed for parallelism ]
!           ========================
            phiexp(n,isph) = phi(ii)
!            
          endif
        enddo
      enddo
! 
!     loop over atoms
      do isph = 1,nsph
!
!       accumulate f -= sum_n U_n^i' Phi_n^i xi(i,n) 
        call fdoga( isph, xi, phiexp, f(:,isph) ) 
!        
      enddo
!
!     initialize index
      i=0
!
!     loop over atoms
      do isph = 1,nsph
!
!       loop over gridpoints
        do n = 1,ngrid
!
!         non-null contribution from grid point
          if ( ui(n,isph).gt.zero ) then
!   
!           advance index
            i=i+1
!
!           compute zeta(i,n)
            zeta(i) = xi(n,isph)*w(n)*ui(n,isph)
!            
          endif
        enddo
      enddo
!
!     time computation of forces
      call system_clock( count = c2 )
!
!     printing
      if ( iprint.gt.0 ) then
!              
        write(*,1010) dble(c2-c1)/dble(cr)
 1010   format(' computation time of ddPCM forces = ',f8.3,' secs.')
! 
      endif
!
!
! =========================  M O D I F Y    H E R E  =========================
!
!     electric field produced by the charges, at the cavity points [ TARGETS ]
      call efld( nsph, charge, csph, ncav, ccav, ef )
!
!     initialize index
      i=0
!
!     loop over atoms
      do isph = 1, nsph
!
!       loop over gridpoints
        do n = 1, ngrid
!
!         non-null contribution from integration point
          if ( ui(n,isph).gt.zero ) then 
!
!           advance index
            i=i+1
!
!           accumulate FIRST contribution to < z , Phi' >
            f(:,isph) = f(:,isph) + zeta(i)*ef(:,i)
!            
          endif
        enddo
      enddo
!
!     electric field produced by the cavity, at the nuclei [ TARGETS ]
      call efld( ncav, zeta, ccav, nsph, csph, ef )
!
!     loop over atoms
      do isph = 1, nsph
!      
!       accumulate SECOND contribution to < z , Phi' >
        f(:,isph) = f(:,isph) + ef(:,isph)*charge(isph)
!        
      enddo
!
! ==========================  E N D    M O D I F Y  ==========================
!
!
!     scale the forces the cosmo factor
      f = 0.5d0*(eps-1.d0)/eps * f
!
!     energy
      e0 = 0.5d0 * sprod( nbasis*nsph, sigma, psi )
!
!     printing
      if ( iprint.ge.2 ) then
!              
        write(*,*)'----------------------------------------------'
        write(*,*)'ddPCM forces (atomic units):'
        write(*,*)''
        write(*,1004)
 1004   format(' atom',13x,'x',13x,'y',13x,'z' )
        do isph = 1,nsph
!
          write(*,1005) isph, f(:,isph)
 1005     format( 1x,i4,3(2x,e12.5) )       
! 
        enddo
        write(*,*)''
        write(*,1006) e0*tokcal
 1006   format(' energy = ',e12.5)     
        write(*,*)'----------------------------------------------'
        write(*,*)''
!        
      endif
!
!
endsubroutine compute_forces
!------------------------------------------------------------------------------      
!
!
!
!------------------------------------------------------------------------------      
! Adjoint problem :
!
!        L^* y = Psi
!   A_eps ^* s = y
!
! PCM problem :
!
!   A_eps G = A_oo F
!       L X = G
!
! Contract :
!
!  < Psi , X' > = < s , A'(F-G) > + 2pi (1 - f_eps) < s , F' > + COSMO
!
! Forces :
!
!   f = -1/2 < Psi , X' >
!
!------------------------------------------------------------------------------      
!
subroutine forces_pcm( Phi, Charge, Psi, Sigma, Flm, F )
!
      use ddcosmo , only : zero, pi, ngrid, nsph, nbasis, lmax, intrhs, basis, &
                           iprint, ncav, fdoga, eps, wghpot, sprod, ui, w, ccav, &
                           csph
!      
      implicit none
      real*8, dimension(ncav),        intent(in)  :: Phi
      real*8, dimension(       nsph), intent(in)  :: Charge
      real*8, dimension(nbasis,nsph), intent(in)  :: Psi
      real*8, dimension(nbasis,nsph), intent(in)  :: Sigma
      real*8, dimension(nbasis,nsph), intent(in)  :: Flm
      real*8, dimension(3,nsph),      intent(out) :: F
!
      real*8, dimension(ngrid,nsph) :: phiexp, xi
      real*8, dimension(nbasis,nsph) :: glm, s, y
      real*8 :: rvoid, e0, xx(1), f_eps, faux(3)
      real*8 :: zeta(ncav), ef(3,ncav)
!
      integer :: isph, n, c1, c2, cr, i
      real*8, parameter :: tokcal=627.509469d0
!
!-------------------------------------------------------------------------------------
!
!     COSMO factor
      f_eps = 0.5d0*(eps-1.d0)/eps
!      
!     initialize the timer
      call system_clock( count_rate = cr )
      call system_clock( count = c1 )
!
!
!     0. Adjoint solves
!     -----------------
!   
!     L^* y = Psi     
      call cosmo( .true., .false., xx, xx, Psi, y, rvoid )
!
!     A_eps^* s = y
      call pcm( .true., .false., .true., xx, y, s )
!
!
!     1. COSMO contribution
!     ---------------------
!
      F = zero
      call forces( nsph, Charge, Phi, Sigma, y, F )
!
!     remove COSMO factor, flip sign
      F = -F/f_eps
!
!
!     2. < , > = < , > + < s , A'(F-G) >
!     ----------------------------------
!
!     expand Phi
      call wghpot( -Phi, phiexp )
!      
!     compute glm
      do isph = 1,nsph
!      
        call intrhs( isph, phiexp(:,isph), glm(:,isph) )
!        
      enddo
!      
!     contract
      do isph = 1,nsph 
!      
        call contract_dA( s, (Flm-glm) , isph, faux(:) )
        F(:,isph) = F(:,isph) + faux(:)
!
      enddo
!
!     3. < , > = < , > - 4pi/(eps-1) < s , F' >
!     ------------------------------------------------
!
!     compute xi
!     ----------
!
!     expand s, rescale, flip sign
      do isph = 1, nsph
        do n = 1, ngrid
!        
          xi(n,isph) = 4.d0*pi/(eps-1.d0) * dot_product( s(:,isph), basis(:,n) )
!          
        enddo
      enddo
!
!     < xi , F' >
      do isph = 1,nsph 
!
        call fdoga( isph, xi, phiexp, F(:,isph) ) 
!      
      enddo
!
!     compute zeta
!     ------------
!
!     initialize index
      i=0
!
!     loop over atoms
      do isph = 1,nsph
!
!       loop over gridpoints
        do n = 1,ngrid
!
!         non-null contribution from grid point
          if ( ui(n,isph).gt.zero ) then
!   
!           advance index
            i=i+1
!
!           compute zeta(i,n)
            zeta(i) = xi(n,isph)*w(n)*ui(n,isph)
!            
          endif
        enddo
      enddo
!      
!     time computation of forces
      call system_clock( count = c2 )
!      
!     printing
      if ( iprint.gt.0 ) then
!              
        write(*,1010) dble(c2-c1)/dble(cr)
 1010   format(' computation time of ddPCM forces = ',f8.3,' secs.')
! 
      endif
!
!
! =========================  M O D I F Y    H E R E  =========================
!
!     electric field produced by the charges, at the cavity points [ TARGETS ]
      call efld( nsph, charge, csph, ncav, ccav, ef )
!
!     initialize index
      i=0
!
!     loop over atoms
      do isph = 1, nsph
!
!       loop over gridpoints
        do n = 1, ngrid
!
!         non-null contribution from integration point
          if ( ui(n,isph).gt.zero ) then 
!
!           advance index
            i=i+1
!
!           accumulate FIRST contribution to < z , Phi' >
            F(:,isph) = F(:,isph) + zeta(i)*ef(:,i)
!            
          endif
        enddo
      enddo
!
!     electric field produced by the cavity, at the nuclei [ TARGETS ]
      call efld( ncav, zeta, ccav, nsph, csph, ef )
!
!     loop over atoms
      do isph = 1, nsph
!      
!       accumulate SECOND contribution to < z , Phi' >
        F(:,isph) = F(:,isph) + ef(:,isph)*charge(isph)
!        
      enddo
!
! ==========================  E N D    M O D I F Y  ==========================
!
!
!     4. Scale the forces the PCM factor, and flip sign
!     -------------------------------------------------
      f_eps = 0.5d0*(eps-1.d0)/eps
      F = -f_eps * F
!
!     energy
      e0 = 0.5d0 * sprod( nbasis*nsph, Sigma, Psi )
!
!     printing
      if ( iprint.ge.2 ) then
!              
        write(*,*)'----------------------------------------------'
        write(*,*)'ddPCM forces (atomic units):'
        write(*,*)''
        write(*,1004)
 1004   format(' atom',13x,'x',13x,'y',13x,'z' )
        do isph = 1,nsph
!
          write(*,1005) isph, f(:,isph)
 1005     format( 1x,i4,3(2x,e12.5) )       
! 
        enddo
        write(*,*)''
        write(*,1006) e0*tokcal
 1006   format(' energy = ',e12.5)     
        write(*,*)'----------------------------------------------'
        write(*,*)''
!        
      endif
!
!
endsubroutine forces_pcm
!------------------------------------------------------------------------------      
!
!
!
!
!------------------------------------------------------------------------------      
subroutine service_routine1_new( s, x, isph, f )
!
      use ddcosmo , only : ui, nsph, nbasis, zero, ngrid, w, one, basis, csph, &
                           rsph, grid, zi, lmax, dbasis, du, pi, two, intrhs
!
      implicit none 
      real*8, dimension(nbasis,nsph), intent(in)    :: s
      real*8, dimension(nbasis,nsph), intent(in)    :: x
      integer,                        intent(in)    :: isph
      real*8, dimension(3),           intent(inout) :: f
!
      real*8,  dimension(nbasis) :: basloc, vplm
      real*8,  dimension(3,nbasis) :: dbsloc
      real*8,  dimension(lmax+1) :: vcos, vsin
      real*8, dimension(3) :: f4,vij,s_ijn,dt_ijn,f2
      real*8, dimension(3,3) :: ds_ijn
      real*8 :: vvij,t_ijn,f3!!,err
      integer :: n,icomp,jcomp,jsph,ksph,ivoid,ind,l,m
      real*8 :: fgrid(3,ngrid), flm(3,nbasis),fac,ss
!      
!------------------------------------------------------------------------------      
!
!     initialize
      f(:) = zero
!
!
!     diagonal term
!     =============
!
!     1st loop over spheres
      do jsph = 1,nsph
!
!       loop over integration points
        do n = 1,ngrid
!
!         initialize
          ss = zero
!
!         contract over l'
          do l = 0,lmax
!
!           build factor : 2pi / (2l+1)
            fac = two*pi/(two*dble(l)+one)
!
!           compute 1st index
            ind = l*l + l + 1
!
!           contract over m'
            do m = -l,l
!
              ss = ss + fac * basis(ind+m,n) * x(ind+m,jsph)          
!
            enddo
          enddo
!
!         multiply by d_i U_j^n
          fgrid(:,n) = du(:,isph,n,jsph)*ss
!
        enddo
!
!       loop over components
        do icomp = 1,3
!
!         integrate against SH
          call intrhs( ivoid, fgrid(icomp,:), flm(icomp,:) )
!
!         contract over l,m
          f(icomp) = f(icomp) + dot_product( flm(icomp,:), s(:,jsph) )

        enddo
      enddo
!
!
!     non-diagonal terms
!     ==================
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize
        f4 = zero  
!      
!       1st loop over spheres
        do jsph = 1,nsph
!
!         contract over l,m
          f3 = dot_product( basis(:,n), s(:,jsph) )
!              
!         2nd loop over spheres
          do ksph = 1,nsph
!
!
!           skip diagonal term
            if ( ksph.eq.jsph )  cycle
!
!
!           compute s_ijn, t_ijn
            vij(:)   = csph(:,jsph) + rsph(jsph)*grid(:,n) - csph(:,ksph)
            vvij     = sqrt( dot_product( vij(:), vij(:) ) )
            t_ijn    = rsph(ksph)/vvij 
            s_ijn(:) =     vij(:)/vvij

!!!            dt_ijn = zero ; ds_ijn = zero
!!!
!!!!           compute derivatives \grad_i t_ijn
!!!            if ( isph.eq.ksph ) then
!!!              dt_ijn(1:3) = rsph(ksph) * vij(1:3) / vvij**3
!!!!              
!!!!             compute derivatives ds_kjn,icomp / dr_j,jcomp = ds_ijn(icomp,jcomp)
!!!              do icomp = 1,3
!!!                do jcomp = 1,3
!!!!
!!!                  ds_ijn(icomp,jcomp) = vij(icomp)*vij(jcomp) / (vvij**3)
!!!!
!!!                  if ( icomp.eq.jcomp ) then
!!!!
!!!                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) - one / vvij
!!!!
!!!                  endif
!!!!                  
!!!                enddo
!!!              enddo
!!!             endif
!!!
!!!            if ( isph.eq.jsph ) then
!!!              dt_ijn(1:3) = -rsph(jsph) * vij(1:3) / vvij**3
!!!!              
!!!!             compute derivatives ds_kjn,icomp / dr_j,jcomp = ds_ijn(icomp,jcomp)
!!!              do icomp = 1,3
!!!                do jcomp = 1,3
!!!!
!!!                  ds_ijn(icomp,jcomp) = -vij(icomp)*vij(jcomp) / (vvij**3)
!!!!
!!!                  if ( icomp.eq.jcomp ) then
!!!!
!!!                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) +  one / vvij
!!!!
!!!                  endif
!!!!                  
!!!                enddo
!!!              enddo
!!!             endif
!
!
!
!           derivatives of s, t
            if ( isph.eq.jsph ) then
!
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = -rsph(ksph) * vij(1:3) / vvij**3
!
!             compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
              do icomp = 1,3
                do jcomp = 1,3
!
                  ds_ijn(icomp,jcomp) = -vij(icomp)*vij(jcomp) / vvij**3
!
                  if ( icomp.eq.jcomp ) then
!
                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) + one / vvij
!
                  endif
!                  
                enddo
              enddo
!              
            elseif ( isph.eq.ksph ) then
!
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = rsph(ksph) * vij(1:3) / vvij**3
!              
!             compute derivatives ds_kjn,icomp / dr_j,jcomp = ds_ijn(icomp,jcomp)
              do icomp = 1,3
                do jcomp = 1,3
!
                  ds_ijn(icomp,jcomp) = vij(icomp)*vij(jcomp) / vvij**3
!
                  if ( icomp.eq.jcomp ) then
!
                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) - one / vvij
!
                  endif
!                  
                enddo
              enddo
!
            else
!                    
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = zero
!              
!             compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
              ds_ijn(1:3,1:3) = zero
!              
            endif
!             
!           compute Y_l'^m'(s_ijn) , d_i Y_l'^m' ( ... )
            call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!           compute f2(j,n)
!           ---------------
            call compute_grad( isph, jsph, n, t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,ksph), f2(1:3) )
!                    
!           accumulate for f4(n)
!           --------------------
            f4(1:3) = f4(1:3) - f3 * f2(1:3)

          enddo
!
        enddo
!
!       accumulate over n
        f(1:3) = f(1:3) + w(n) * f4(1:3)
!
      enddo
!                      
!
endsubroutine service_routine1_new
!----------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------
subroutine compute_grad( isph, jsph, n, t, dt, ds, basloc, dbsloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, four, pi, ui, du, zi
!
      implicit none
      integer,                     intent(in)  :: isph
      integer,                     intent(in)  :: jsph
      integer,                     intent(in)  :: n
      real*8,                      intent(in)  :: t
      real*8, dimension(3  ),      intent(in)  :: dt
      real*8, dimension(3,3),      intent(in)  :: ds
      real*8, dimension(  nbasis), intent(in)  :: basloc
      real*8, dimension(3,nbasis), intent(in)  :: dbsloc
      real*8, dimension(  nbasis), intent(in)  :: x
      real*8, dimension(3),        intent(out) :: f2
!
      real*8, dimension(3,nbasis) :: s1,s2,f1
      real*8, dimension(3) :: s3
      real*8 :: tt,fl,fac
      integer :: icomp,jcomp,l,m,ind
!      
!----------------------------------------------------------------------------------
! Recall :
!
!   f1(1:3,lm) = \grad [ U * t^l+1 * Y(s) ]
!
!              = (\grad U) t^l+1 * Y(s) + U * t^l [ (l+1) (\grad t) Y(s) + t ( \grad s )^T \grad Y(s) ]
!
!              = t^l [ (l+1) (\grad t) U * Y + t [ (\grad U) Y + U ( \grad s )^T \grad Y ] ]
!
!              = tt  [ (l+1) s1(1:3,lm)      + t   s2(1:3,lm) ]
!
! and :
!
!                  4 pi l
!   f2(1:3) = sum  ------  sum x(lm) * f1(1:3,lm)
!              l   2l + 1   m
!
!                  4 pi l
!           = sum  ------  s3(1:3,l)  
!              l   2l + 1 
!
!           = sum fac(l) * s3(1:3,l) 
!              l
!
!----------------------------------------------------------------------------------
!
!     initialize f2
      f2(1:3) = zero
!
!     compute s1
      do icomp = 1,3
!      
        s1(icomp,1:nbasis) = dt(icomp) * basloc(1:nbasis) * ui(n,jsph)
!        
      enddo
!
!     compute s2
      s2(1:3,1:nbasis) = zero
!      
      do icomp = 1,3
!
!       accumulate
        do jcomp = 1,3
!
          s2(icomp,1:nbasis) = s2(icomp,1:nbasis) + dbsloc(jcomp,1:nbasis)*ds(jcomp,icomp)
!
        enddo
!
        s2(icomp,1:nbasis) = ui(n,jsph)*s2(icomp,1:nbasis) + du(icomp,isph,n,jsph)*basloc(1:nbasis)
!
      enddo
!
!     initialize factor t^l
      tt = one
!
!     contract over l
      do l = 0,lmax
!
!       build factor : 4pi*l / (2l+1)
        fl = dble(l)
        fac = four*pi*fl / (two*fl+one)
!
!       build f1
        f1(1:3,1:nbasis) = tt * ( (fl+one)*s1(1:3,1:nbasis) + t*s2(1:3,1:nbasis) )
!
!       compute 1st index
        ind = l*l + l + 1
!
!       initialize s3(l)
        s3(1:3) = zero
!
!       contract over m
        do m = -l,l
!
          s3(1:3) = s3(1:3) + f1(1:3,ind+m) * x(ind+m)          
!
        enddo
!
!       accumulate for f2
        f2(1:3) = f2(1:3) + fac*s3(1:3)
!
!       update factor t^l
        tt = tt*t
!
      enddo
!
!
endsubroutine compute_grad
!------------------------------------------------------------------------------      
!
!
!
!------------------------------------------------------------------------------      
subroutine contract_dA( s, x, isph, f )
!
      use ddcosmo , only : ui, nsph, nbasis, zero, ngrid, w, one, basis, csph, &
                           rsph, grid, zi, lmax, dbasis, du, pi, two, intrhs, &
                           inl, nl
!
      implicit none 
      real*8, dimension(nbasis,nsph), intent(in)    :: s
      real*8, dimension(nbasis,nsph), intent(in)    :: x
      integer,                        intent(in)    :: isph
      real*8, dimension(3),           intent(inout) :: f
!
      real*8,  dimension(nbasis) :: basloc, vplm
      real*8,  dimension(3,nbasis) :: dbsloc
      real*8,  dimension(lmax+1) :: vcos, vsin
      real*8, dimension(3) :: f4,vij,s_ijn,dt_ijn,f2
      real*8, dimension(3,3) :: ds_ijn
      real*8 :: vvij,t_ijn,f3!!,err
      integer :: n,icomp,jcomp,jsph,ksph,ivoid,ind,l,m,i
      real*8 :: fgrid(3,ngrid), flm(3,nbasis),fac,ss
!      
!------------------------------------------------------------------------------      
!
!
!!!!     s_i ( sum grad_i A_ik * x_k )
!!!!            k
!!!!
!!!      do n = 1,ngrid
!!!!
!!!        w(n)*basis(:,n)
!!!
!!!        do ksph = 1,nsph
!!!
!!!        enddo
!!!      enddo
!!!!
!!!!
!!!!     x_i (   sum   s_j * grad_i A_ji )
!!!!           j \ne i
!!!!
!!!      do jsph = 1,nsph
!!!!
!!!        if ( jsph.eq.isph )  cycle
!!!
!!!      enddo
!!!
!!!!       sum   x_k (    sum    s_j * grad_i A_jk )
!!!!     k \ne i       j \in N_i
!!!!
!!!      do ksph = 1,nsph
!!!!      
!!!        if ( ksph.eq.isph )  cycle
!!!
!!!        do i = inl(isph),inl(isph+1)-1
!!!
!!!          jsph = nl(i)
!!!
!!!        enddo
!!!      enddo

!
      f=zero
!
!     diagonal term
!     =============
!
!     1st loop over spheres
      do jsph = 1,nsph
!
!       loop over integration points
        do n = 1,ngrid
!
!         initialize
          ss = zero
!
!         contract over l'
          do l = 0,lmax
!
!           build factor : 2pi / (2l+1)
            fac = two*pi/(two*dble(l)+one)
!
!           compute 1st index
            ind = l*l + l + 1
!
!           contract over m'
            do m = -l,l
!
              ss = ss + fac * basis(ind+m,n) * x(ind+m,jsph)          
!
            enddo
          enddo
!
!         multiply by d_i U_j^n
          fgrid(:,n) = du(:,isph,n,jsph)*ss
!
        enddo
!
!       loop over components
        do icomp = 1,3
!
!         integrate against SH
          call intrhs( ivoid, fgrid(icomp,:), flm(icomp,:) )
!
!         contract over l,m
          f(icomp) = f(icomp) + dot_product( flm(icomp,:), s(:,jsph) )

        enddo
      enddo
!
!
!     non-diagonal terms
!     ==================
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize
        f4 = zero  
!              
!       1nd loop over spheres
        do ksph = 1,nsph
!
!
!         2nd loop over (proper!) neighbors
          do i = inl(isph),inl(isph+1)-1
!
            jsph = nl(i)
!
!           skip diagonal term
            if ( ksph.eq.jsph )  cycle
!
!           contract over l,m
            f3 = dot_product( basis(:,n), s(:,jsph) )
!
!           compute s_ijn, t_ijn
            vij(:)   = csph(:,jsph) + rsph(jsph)*grid(:,n) - csph(:,ksph)
            vvij     = sqrt( dot_product( vij(:), vij(:) ) )
            t_ijn    = rsph(ksph)/vvij 
            s_ijn(:) =     vij(:)/vvij
!
!           derivatives of s, t
            if ( isph.eq.jsph ) then
!
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = -rsph(ksph) * vij(1:3) / vvij**3
!
!             compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
              do icomp = 1,3
                do jcomp = 1,3
!
                  ds_ijn(icomp,jcomp) = -vij(icomp)*vij(jcomp) / vvij**3
!
                  if ( icomp.eq.jcomp ) then
!
                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) + one / vvij
!
                  endif
!                  
                enddo
              enddo
!              
            elseif ( isph.eq.ksph ) then
!
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = rsph(ksph) * vij(1:3) / vvij**3
!              
!             compute derivatives ds_kjn,icomp / dr_j,jcomp = ds_ijn(icomp,jcomp)
              do icomp = 1,3
                do jcomp = 1,3
!
                  ds_ijn(icomp,jcomp) = vij(icomp)*vij(jcomp) / vvij**3
!
                  if ( icomp.eq.jcomp ) then
!
                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) - one / vvij
!
                  endif
!                  
                enddo
              enddo
!
            else
!                    
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = zero
!              
!             compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
              ds_ijn(1:3,1:3) = zero
!              
            endif
!             
!           compute Y_l'^m'(s_ijn) , d_i Y_l'^m' ( ... )
            call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!           compute f2(j,n)
!           ---------------
            call compute_grad( isph, jsph, n, t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,ksph), f2(1:3) )
!                    
!           accumulate for f4(n)
!           --------------------
            f4(1:3) = f4(1:3) - f3 * f2(1:3)

          enddo
!
!
!         case j = i
          jsph = isph
!
!         skip diagonal term
          if ( ksph.eq.jsph )  cycle
!
!         contract over l,m
          f3 = dot_product( basis(:,n), s(:,jsph) )
!
!         compute s_ijn, t_ijn
          vij(:)   = csph(:,jsph) + rsph(jsph)*grid(:,n) - csph(:,ksph)
          vvij     = sqrt( dot_product( vij(:), vij(:) ) )
          t_ijn    = rsph(ksph)/vvij 
          s_ijn(:) =     vij(:)/vvij
!
!         derivatives of s, t
          if ( isph.eq.jsph ) then
!
!           compute derivatives \grad_i t_ijn
            dt_ijn(1:3) = -rsph(ksph) * vij(1:3) / vvij**3
!
!           compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
            do icomp = 1,3
              do jcomp = 1,3
!
                ds_ijn(icomp,jcomp) = -vij(icomp)*vij(jcomp) / vvij**3
!
                if ( icomp.eq.jcomp ) then
!
                  ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) + one / vvij
!
                endif
!                
              enddo
            enddo
!            
          elseif ( isph.eq.ksph ) then
!
!           compute derivatives \grad_i t_ijn
            dt_ijn(1:3) = rsph(ksph) * vij(1:3) / vvij**3
!            
!           compute derivatives ds_kjn,icomp / dr_j,jcomp = ds_ijn(icomp,jcomp)
            do icomp = 1,3
              do jcomp = 1,3
!
                ds_ijn(icomp,jcomp) = vij(icomp)*vij(jcomp) / vvij**3
!
                if ( icomp.eq.jcomp ) then
!
                  ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) - one / vvij
!
                endif
!                
              enddo
            enddo
!
          else
!                  
!           compute derivatives \grad_i t_ijn
            dt_ijn(1:3) = zero
!            
!           compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
            ds_ijn(1:3,1:3) = zero
!            
          endif
!           
!         compute Y_l'^m'(s_ijn) , d_i Y_l'^m' ( ... )
          call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!         compute f2(j,n)
!         ---------------
          call compute_grad( isph, jsph, n, t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,ksph), f2(1:3) )
!                  
!         accumulate for f4(n)
!         --------------------
          f4(1:3) = f4(1:3) - f3 * f2(1:3)



        enddo
!
!       accumulate over n
        f(1:3) = f(1:3) + w(n) * f4(1:3)
!
      enddo
!                      
!
endsubroutine contract_dA
