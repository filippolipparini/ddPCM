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
      use ddcosmo , only : zero, ngrid, nsph, nbasis, zero, lmax, intrhs, itsolv2, &
                           basis, fdoka, fdokb, iprint, ncav, ui, ccav, csph, &
                           w, fdoga, eps, sprod
!      
      implicit none
      real*8, dimension( ngrid,nsph), intent(in)  :: Phi
      real*8, dimension(       nsph), intent(in)  :: charge
      real*8, dimension(nbasis,nsph), intent(in)  :: Psi
      real*8, dimension(nbasis,nsph), intent(in)  :: sigma
      real*8, dimension(nbasis,nsph), intent(in)  :: Phi_eps
      real*8, dimension(nsph,3),      intent(out) :: f
!
      real*8, dimension(ngrid,nsph) :: xi
      real*8, dimension(nbasis,nsph) :: w_lm, s, y, z
      real*8, dimension(ngrid) :: x
      real*8, dimension(nbasis) :: xlm, basloc, vplm
      real*8, dimension(3,nbasis) :: dbsloc
      real*8, dimension(lmax+1) :: vcos, vsin
      real*8 :: rvoid,e0
      real*8 :: ef(3,ncav),zeta(ncav)
!
      integer :: isph, jsph, icomp, n, i, c1, c2, cr
      logical, parameter :: star=.true.
      real*8, parameter :: tokcal=627.509469d0
!
!-------------------------------------------------------------------------------------
!
!     initialize
      f(:,:) = zero
!
!     initialize the timer
      call system_clock( count_rate = cr )
      call system_clock( count = c1 )
!      
!
!     STEP 1 : solve adjoint problem (A_eps L)^T s = Psi
!     --------------------------------------------------
!   
!     solve L^T y = Psi     
      call itsolv2( star, .true., Psi, Psi, y, rvoid )
!
!     solve A_eps^T s = y
      call ADJpcm( y, s )
!
!
!     STEP 2 : compute f = - < s , A' ( Phi - Phi_eps ) >
!     ---------------------------------------------------
!
!     compute SH expansion of Phi on i-sphere
      do isph = 1,nsph
!      
        call intrhs( isph, Phi(:,isph), w_lm(:,isph) )
!        
      enddo
!
!     compute w = Phi - Phi_eps
      w_lm(:,:) = w_lm(:,:) - Phi_eps(:,:)
!
!     contract
      do isph = 1,nsph 
!      
        call service_routine1( s, w_lm , isph, f(isph,1:3) )
!
      enddo
!
!     flip sign 
      f = -f
!
!
      write(*,*)'----------------------------------------------'
      write(*,*)'ddPCM forces (atomic units): [1]'
      write(*,*)''
      write(*,1004)
      do isph = 1,nsph
!
        write(*,1005) isph, f(isph,:)
! 
      enddo
      write(*,*)'----------------------------------------------'
      write(*,*)''
!!!      read(*,*)
!
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
!     initialize
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
!     loop over atoms
      do isph = 1,nsph
!
!       accumulate f -= sum_n U_n^i' Phi_n^i xi(i,n) 
        call fdoga( isph, xi, Phi, f(isph,:) ) 
!        
      enddo
!
!
      write(*,*)'----------------------------------------------'
      write(*,*)'ddPCM forces (atomic units): [2]'
      write(*,*)''
      write(*,1004)
      do isph = 1,nsph
!
        write(*,1005) isph, f(isph,:)
! 
      enddo
      write(*,*)'----------------------------------------------'
      write(*,*)''
!!!      read(*,*)
!
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
            f(isph,:) = f(isph,:) + zeta(i)*ef(:,i)
!            
          endif
        enddo
      enddo
!
!
      write(*,*)'----------------------------------------------'
      write(*,*)'ddPCM forces (atomic units): [3]'
      write(*,*)''
      write(*,1004)
      do isph = 1,nsph
!
        write(*,1005) isph, f(isph,:)
! 
      enddo
      write(*,*)'----------------------------------------------'
      write(*,*)''
!!!      read(*,*)
!
!
!     electric field produced by the cavity, at the nuclei [ TARGETS ]
      call efld( ncav, zeta, ccav, nsph, csph, ef )
!
!     loop over atoms
      do isph = 1, nsph
!      
!       accumulate SECOND contribution to < z , Phi' >
        f(isph,:) = f(isph,:) + ef(:,isph)*charge(isph)
!        
      enddo
!
! ==========================  E N D    M O D I F Y  ==========================
!
!
      write(*,*)'----------------------------------------------'
      write(*,*)'ddPCM forces (atomic units): [4]'
      write(*,*)''
      write(*,1004)
      do isph = 1,nsph
!
        write(*,1005) isph, f(isph,:)
! 
      enddo
      write(*,*)'----------------------------------------------'
      write(*,*)''
!!!      read(*,*)
!
!
!
!     STEP 4 : compute f = f + < y , L' sigma >
!     -----------------------------------------
!
!     initialize
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
        call fdoka( isph, sigma, xi(:,isph), basloc, dbsloc, vplm, vcos, vsin, f(isph,:) ) 
        call fdokb( isph, sigma, xi,         basloc, dbsloc, vplm, vcos, vsin, f(isph,:) ) 
!
      enddo
!
!     scale the forces the cosmo factor
      f = 0.5d0*(eps-1.d0)/eps * f
!
!
!     time computation of forces
      call system_clock( count = c2 )
!
!     printing
      if ( iprint.gt.0 ) then
!              
        write(*,1010) dble(c2-c1)/dble(cr)
 1010   format(' computation of ddPCM forces : ',f8.3,' sec')
! 
      endif
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
          write(*,1005) isph, f(isph,:)
 1005     format( 1x,i4,3(2x,e12.5) )       
! 
        enddo
        write(*,*)''
        write(*,1006) e0*tokcal
 1006   format(' energy = ',e12.5)     
        write(*,*)'----------------------------------------------'
        write(*,*)''
!!!        read(*,*)
!        
      endif
!
!
endsubroutine compute_forces
!------------------------------------------------------------------------------      
!
!
!
!
!------------------------------------------------------------------------------      
! Compute contraction :
!
!   f =+ < s , dA/dr_i x > = 
!
!     =+ sum s_j ( sum (d_i A)_jk x_k )
!         j         k
!
! Recall that (d_i A)_jj = 0 for i \ne j . Then :
!
!     =+ sum s_j (   sum  (d_i A)_jk x_k + (d_i A)_jj x_j )
!         j        k \ne j
!
!     =+ sum s_j    sum  (d_i A)_jk x_k + sum s_j (d_i A)_jj x_j
!         j       k \ne j                  j
!
!     =+ sum s_j    sum  (d_i A)_jk x_k + s_i (d_i A)_ii x_i
!         j       k \ne j 
!
!     =+   sum   s_j   sum   (d_i A)_jk x_k + s_i   sum   (d_i A)_ik x_k +  s_i (d_i A)_ii x_i
!        j \ne i     k \ne j                      k \ne i
!
! Finally, since (d_i A)_jk = 0 for k \ne i and j \ne i , then :
!
!     =+   sum   s_j (d_i A)_ji x_i + s_i   sum   (d_i A)_ik x_k +  s_i (d_i A)_ii x_i
!        j \ne i                          k \ne i
!
!             +---+                       +--------------------+        +---+
!             |///|                       |///////////|  |/////|        |///|
!             |///|                       +--------------------+        +---+
!             |///|
!             |///|
!             |///|
!             |///|
!             |---|
!             |   |
!             |---|
!             |///|
!             +---+
!------------------------------------------------------------------------------      
subroutine service_routine1( s, x, isph, f )
!
      use ddcosmo , only : ui, nsph, nbasis, zero, ngrid, w, one, basis, csph, &
                           rsph, grid, zi, lmax, dbasis
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
      integer :: n,icomp,jcomp,jsph,ksph
!      
!------------------------------------------------------------------------------      
!
!     case j \ne i : I_1
!     ==================
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize f4(n)
        f4(1:3) = zero
!
!       1st loop over spheres
        do jsph = 1,nsph
!
          if ( jsph.ne.isph ) then
!
!           non-zero contribution
            if ( ui(n,jsph).gt.zero ) then
!
!             compute f3(j,n)
!             ---------------
              f3 = ui(n,jsph) * dot_product( basis(:,n), s(:,jsph) )
!
!             compute s_ijn, t_ijn
              vij(:)   = csph(:,jsph) + rsph(jsph)*grid(:,n) - csph(:,isph)
              vvij     = sqrt( dot_product( vij(:), vij(:) ) )
              t_ijn    = rsph(isph)/vvij 
              s_ijn(:) =     vij(:)/vvij
!
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = rsph(isph) * vij(1:3) / vvij**3
!
!             compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
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

!!!              err=zero
!!!              do icomp=1,3
!!!              do jcomp=1,3
!!!                err = err + (ds_ijn(icomp,jcomp)-ds_ijn(jcomp,icomp))**2
!!!              enddo
!!!              enddo
!!!              write(*,*)'err = ',sqrt(err)

!             
!             compute Y_ll^mm(s_ijn)
              call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!             compute f2(j,n)
!             ---------------
              call compute_f2_jn( t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,isph), f2(1:3) )
!                      
!             accumulate for f4(n)
!             --------------------
              f4(1:3) = f4(1:3) + f2(1:3) * f3
!                      
            endif  
          endif  
        enddo  
!
!       accumulate for f
!       ----------------
        f(1:3) = f(1:3) - w(n) * f4(1:3)
!
      enddo  
!
!
!     case j \eq i : I_2,1
!     ====================
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize f4(n)
        f4(1:3) = zero
!
!       1st loop over spheres
        do ksph = 1,nsph
!
          if ( ksph.ne.isph ) then
!
!           compute f3(k,n)
!           ---------------
            f3 = dot_product( basis(:,n), s(:,isph) )
!
!           compute s_kjn, t_kjn
            vij    = csph(:,isph) + rsph(isph)*grid(:,n) - csph(:,ksph)
            vvij   = sqrt(dot_product(vij,vij))
            t_ijn  = rsph(ksph)/vvij 
            s_ijn  =        vij/vvij
!
!           compute derivatives \grad_j t_kjn
            dt_ijn(1:3) = - rsph(ksph) * vij(1:3) / vvij**3
!
!           compute derivatives ds_kjn,icomp / dr_j,jcomp = ds_ijn(icomp,jcomp)
            do icomp = 1,3
              do jcomp = 1,3
!
                ds_ijn(icomp,jcomp) = - vij(icomp)*vij(jcomp) / vvij**3
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
!           compute Y_ll^mm(s_kjn)
            call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!           compute f2(j,n)
!           ---------------
            call compute_f2_kn( isph, n, t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,ksph), f2(1:3) )
!                    
!           accumulate for f4(n)
!           --------------------
            f4(1:3) = f4(1:3) + f2(1:3) * f3
!                      
          endif  
        enddo  
!
!       accumulate for f
!       ----------------
        f(1:3) = f(1:3) - w(n) * f4(1:3)
!
      enddo  
!
!
!     case j \eq i : I_2,2
!     ====================
!
!     loop over integration points
      do n = 1,ngrid
!
!       compute f3(n) [ store in f2 for convenience ]
!       ---------------------------------------------
        f2(1:3) = zi(1:3,n,isph) * dot_product( basis(:,n), s(:,isph) )
!
!       compute f2(n) [ store in f3 for convenience ]
!       ---------------------------------------------
        call compute_f2_n( basloc, x(:,isph), f3 )
!
!       accumulate for f
!       ----------------
        f(1:3) = f(1:3) + w(n) * f2(1:3) * f3
!
      enddo  
!
!
endsubroutine service_routine1
!----------------------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------------------
subroutine compute_f2_jn( t, dt, ds, basloc, dbsloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, four, pi
!
      implicit none
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
!   f1(1:3,lm) = \grad [ t^l+1 * Y(s) ]
!
!              = t^l  [ (l+1) \grad t Y(s) + t ( \grad s )^T \grad Y(s) ]
!
!              = tt   [ (l+1) s1(1:3,lm)   + s2(1:3,lm)                 ]
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
        s1(icomp,1:nbasis) = dt(icomp)*basloc(1:nbasis)
!        
      enddo
!
!     compute s2
      do icomp = 1,3
!
        s2(icomp,1:nbasis) = zero
!
!       accumulate
        do jcomp = 1,3
!
          s2(icomp,1:nbasis) = s2(icomp,1:nbasis) + dbsloc(jcomp,1:nbasis)*ds(jcomp,icomp)
!
        enddo
!
        s2(icomp,1:nbasis) = t * s2(icomp,1:nbasis) 
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
        f1(1:3,1:nbasis) = tt * ( (fl+one)*s1(1:3,1:nbasis) + s2(1:3,1:nbasis) )
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
endsubroutine compute_f2_jn
!----------------------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------------------
subroutine compute_f2_kn( jsph, n, t, dt, ds, basloc, dbsloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, four, pi, ui, zi
!
      implicit none
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
        s2(icomp,1:nbasis) = ui(n,jsph)*s2(icomp,1:nbasis) + zi(icomp,n,jsph)*basloc(1:nbasis)
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
endsubroutine compute_f2_kn
!----------------------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------------------
subroutine compute_f2_n( basloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, pi
!
      implicit none
      real*8, dimension(nbasis), intent(in)  :: basloc
      real*8, dimension(nbasis), intent(in)  :: x
      real*8,                    intent(out) :: f2
!
      real*8 :: fl,s,fac
      integer :: l,m,ind
!      
!----------------------------------------------------------------------------------
! Recall :
!
!             2 pi 
!   f2 = sum  ----  sum Y(lm) * x(lm)
!         l   2l+1   m
!
!             2 pi
!      = sum  ----  s(l)
!         l   2l+1 
!
!----------------------------------------------------------------------------------
!
!     initialize f2
      f2 = zero
!
!     contract over l
      do l = 0,lmax
!
!       build factor : 2pi / (2l+1)
        fl = dble(l)
        fac = two*pi/(two*fl+one)
!
!       compute 1st index
        ind = l*l + l + 1
!
!       initialize s
        s = zero
!        
!       contract over m
        do m = -l,l
!
          s = s + basloc(ind+m) * x(ind+m)          
!
        enddo
!
!       accumulate for f2
        f2 = f2 + fac*s
!
      enddo
!
!
endsubroutine compute_f2_n
!----------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------
subroutine compute_dphi( dphi )
!
      use ddcosmo , only : ngrid, nsph
!
      real*8, dimension(ngrid,nsph,nsph,3), intent(out) :: dphi
      integer, save :: iflag = 0
!
!----------------------------------------------------------------------------------
!
      if ( iflag.eq.0 ) then
        write(*,*)'compute_dphi : DUMMY routine !'      
        write(*,*)''
        iflag=1
      endif

!     nothing for now ...   
      dphi(:,:,:,:) = 0.d0
!
!
endsubroutine compute_dphi
!----------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------
subroutine ADJcheck
!       
      use ddcosmo , only : nbasis,nsph,ngrid,lmax,zero,one,two,csph,rsph,memfree, &
                           ddinit,ui,zi,iquiet,eps
!
      implicit none
      real*8 :: f(nbasis,nsph),Af( nbasis)
      real*8 :: e(nbasis,nsph),ATe(nbasis)
      real*8 :: s(nbasis,nsph)
      real*8 :: A(nbasis*nsph,nbasis*nsph),AT(nbasis*nsph,nbasis*nsph)
      real*8 :: dA1(nbasis*nsph,nbasis*nsph),dA2(nbasis*nsph,nbasis*nsph)
      real*8 :: dA3(nbasis*nsph,nbasis*nsph)
      real*8 :: dA(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: A_plus(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: xlm(nbasis),xx(ngrid),vplm(nbasis),vcos(lmax+1),vsin(lmax+1), &
                basloc(nbasis)
      integer :: isph,jsph,i,j,ibeg,iend,nsph_save,icomp,ksph,iter,n
      real*8 :: s1,s2,eeps,err,rnorm
      integer, parameter :: niter = 6
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), s3(3)
      real*8 :: x(nsph), y(nsph), z(nsph), iwork(nsph*3,2), rwork(niter,nsph*3)
      real*8 :: rrate(niter,nsph*3)
!
!--------------------------------------------------------------------------------
!
!     activate quiet flag
      iquiet = .true.
!
!     check A = ( A^T )^T
!     -------------------
!
!!!     set epsilon
!!!      eps=two
!
!     construct A, A^T
      do isph = 1,nsph
!
        do jsph = 1,nsph
          do j = 1,nbasis
!
!           standard basis vector e_j
            e(:,   :) = zero
            e(j,jsph) = one
!
!           compute A  _i e_j
            ibeg = (isph-1)*nbasis+1
            iend = (isph-1)*nbasis+nbasis
            call mkrvec( isph, eps, e, A( ibeg:iend,(jsph-1)*nbasis+j), xlm, xx, basloc, vplm, vcos, vsin )
!
!           compute A^T_i e_j
            call ADJvec( isph, eps, e, AT(ibeg:iend,(jsph-1)*nbasis+j), xlm, xx, basloc, vplm, vcos, vsin )
!        
          enddo
        enddo
      enddo
!
!     initialize
      s1 = zero ; s2 = zero
!      
!     compute Frobenious norm of A and A^T
      do isph = 1,nsph
        do i = 1,nbasis
          do jsph = 1,nsph
            do j = 1,nbasis
!
!             accumulate a_ij^2, (a^T)_ij^2
              s1 = s1 + A ( (isph-1)*nbasis+i , (jsph-1)*nbasis+j )**2
              s2 = s2 + AT( (isph-1)*nbasis+i , (jsph-1)*nbasis+j )**2
!
            enddo
          enddo
        enddo
      enddo
!
!     || A^T ||_F     || A ||_F
      s1 = sqrt(s1) ; s2 = sqrt(s2)
!      
!     print
      write(*,1002) abs(s1-s2) / abs(s1)
 1002 format(' | ||A^T||_F - ||A||_F | / ||A||_F = ', e12.5)
      write(*,*) ''
!
!!!      write(*,*) 'A = '
!!!      do isph = 1,nsph
!!!        do i = 1,nbasis
!!!          write(*,"(4x,300(e12.5,2x))") ( A((isph-1)*nbasis+i,j), j=1,nbasis*nsph )
!!!        enddo
!!!      enddo
!!!      write(*,*)''
!!!!
!!!      write(*,*) '(A^T)^T = '
!!!      do isph = 1,nsph
!!!        do i = 1,nbasis
!!!          write(*,"(4x,300(e12.5,2x))") ( AT(j,(isph-1)*nbasis+i), j=1,nbasis*nsph )
!!!        enddo
!!!      enddo
!!!      write(*,*)''
!
!
!!!!     check < A^T e , f > = < e , A f >
!!!!     ---------------------------------
!!!!
!!!!     initialize random number generator
!!!      call random_seed
!!!!      
!!!!     build e, f
!!!      s1=zero ; s2=zero
!!!      do isph = 1,nsph
!!!        do j = 1,nbasis
!!!!        
!!!          call random_number( e(j,isph) )
!!!          call random_number( f(j,isph) )
!!!!          
!!!          s1 = s1 + e(j,isph)**2
!!!          s2 = s2 + f(j,isph)**2
!!!!          
!!!        enddo
!!!      enddo
!!!      e(:,:)=e(:,:)/sqrt(s1)
!!!      f(:,:)=f(:,:)/sqrt(s2)
!!!!
!!!!     initialize
!!!      Af( :)=zero
!!!      ATe(:)=zero
!!!      s1=zero
!!!      s2=zero
!!!!      
!!!      do isph = 1,nsph
!!!!
!!!!       compute A_i f 
!!!        call mkrvec( isph, eps, f, Af( :), xlm, xx, basloc, vplm, vcos, vsin )
!!!!
!!!!       compute A^T_i e
!!!        call ADJvec( isph, eps, e, ATe(:), xlm, xx, basloc, vplm, vcos, vsin )
!!!!
!!!!       accumulate < e, A f >
!!!        s1 = s1 + dot_product( e(:,isph), Af( :) )
!!!!        
!!!!       accumulate < A^T e, f >
!!!        s2 = s2 + dot_product( f(:,isph), ATe(:) )
!!!!        
!!!      enddo
!!!!      
!!!      write(*,1000) abs(s1-s2) / abs(s1)
!!! 1000 format(' | <e,Af> - <A^T f,e> | / |<e,Af>| = ', e12.5)
!!!      write(*,*) ''
!
!
!     check solution of adjoint problem : < e , f > = < A^T s , f > = < s , A f >
!     ---------------------------------------------------------------------------
!
!     initialize random number generator
      call random_seed
!      
!     build e, f 
      s1=zero ; s2=zero
      do isph = 1,nsph
        do j = 1,nbasis
!        
          call random_number( e(j,isph) )
          call random_number( f(j,isph) )
!          
          s1 = s1 + e(j,isph)**2
          s2 = s2 + f(j,isph)**2
!          
        enddo
      enddo
!
!     normalize
      e(:,:)=e(:,:)/sqrt(s1)
      f(:,:)=f(:,:)/sqrt(s2)
!
!     solve A_eps^T s = e
      call ADJpcm( e, s )
!      
!     initialize
      Af(:) = zero ; s1 = zero ; s2 = zero
!      
      do isph = 1,nsph
!
!       compute A_i f 
        call mkrvec( isph, eps, f, Af(:), xlm, xx, basloc, vplm, vcos, vsin )
!
!       accumulate < s , A f >
        s1 = s1 + dot_product( s(:,isph), Af(:) )
!        
!       accumulate < e , f >
        s2 = s2 + dot_product( e(:,isph), f(:,isph) )
!        
      enddo
!
!     print
      write(*,1001) abs(s2-s1) / abs(s1)
 1001 format(' | <e,f> - <s,Af> | / |<s,Af>| = ', e12.5)
      write(*,*) ''
!
!     initialize
      rwork = zero ; rrate = zero
!
!     check derivatives
!     -----------------
!
!     1. compute analytical derivatives
      do ksph = 1,nsph 
        do isph = 1,nsph
          do i = 1,nbasis
            do jsph = 1,nsph
              do j = 1,nbasis
!
!               standard basis vectors e_i, f_j
                e(:,:   )=zero
                f(:,:   )=zero
                e(i,isph)=one
                f(j,jsph)=one
!
!               compute < e, dA/dr_k f > 
                s3=zero
                call service_routine1( e(:,:), f(:,:), ksph, s3 )
!
!               store
                do icomp = 1,3
!
                  dA( (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp ) = s3(icomp)
!            
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
!     save initial DS
      nsph_save = nsph
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
!     set initial increment
      eeps=0.1d0
!
!     loop over increments
      do iter = 1,niter
!        
!       compute A^+
        do ksph = 1,nsph
!
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph) + eeps
            case(2) ; y(ksph) = y_save(ksph) + eeps
            case(3) ; z(ksph) = z_save(ksph) + eeps
            endselect
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!
!           initialize
            err=zero
            rnorm=zero
            A_plus(:,:,ksph,icomp)=zero
!
!           build A^+
            do isph = 1,nsph
              do jsph = 1,nsph
                do j = 1,nbasis
!
                  e(:,:   )=zero
                  e(j,jsph)=one
!
!                 compute A^+_i e_j
                  ibeg = (isph-1)*nbasis+1
                  iend = (isph-1)*nbasis+nbasis
                  call mkrvec( isph, eps, e, A_plus( ibeg:iend,(jsph-1)*nbasis+j,ksph,icomp ), xlm, xx, basloc, vplm, vcos, vsin )
!
!!!                  if ( ( isph.eq.2 ) .and. ( jsph.eq.1 ) ) then

!                 accumulate error
                  do i = 1,nbasis
                    err = err + ( ( A_plus((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) -           &
                                    A(     (isph-1)*nbasis+i,(jsph-1)*nbasis+j           ) ) / eeps -  &
                                    dA(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) )**2
                    rnorm = rnorm + (  dA(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) )**2
!
                  enddo
                
!!!          endif

                enddo
              enddo
            enddo
!            
!!!!           numerical derivatives
!!!            write(*,1007) ksph,icomp
!!! 1007       format(' A(r+r_',i2,',',i1,') =')
!!!! 
!!!            do isph = 1,nsph
!!!              do i = 1,nbasis
!!!                write(*,"(4x,300(e12.5,2x))") ( A_plus((isph-1)*nbasis+i,j,ksph,icomp) , j=1,nbasis*nsph )
!!!              enddo
!!!            enddo
!!!            write(*,*)''
!!!!            
!!!            write(*,1008) ksph,icomp
!!! 1008       format('( A(r+r_',i2,',',i1,') - A(r) ) / eps =')
!!!! 
!!!            do isph = 1,nsph
!!!              do i = 1,nbasis
!!!                write(*,"(4x,300(e12.5,2x))") ( ( A_plus((isph-1)*nbasis+i,j,ksph,icomp) - &
!!!                                                     A((isph-1)*nbasis+i,j))/eeps , j=1,nbasis*nsph )
!!!              enddo
!!!            enddo
!!!            write(*,*)''
!
!           store relative error
            rwork(iter,(ksph-1)*3+icomp) = sqrt( err / rnorm )
!
!           store rate of convergence
            if ( iter.gt.1 ) then 
              rrate(iter,(ksph-1)*3+icomp) =  log( rwork(iter-1,(ksph-1)*3+icomp) / &
                                                   rwork(iter  ,(ksph-1)*3+icomp)   ) / log(0.5d0)  
            endif
!
          enddo
        enddo
!
!       update increment
        eeps = eeps/2.d0
!        
      enddo

      eeps = eeps*2.d0
!
!     printing relative error
      write(*,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dA / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''

!     printing rate of convergence
      write(*,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dA / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''

!
!!!!     printing
!!!      do ksph = 1,nsph
!!!        do icomp = 1,3
!!!!        
!!!!         analytical derivatives
!!!          write(*,1005) ksph,icomp
!!! 1005     format(' dA / dr_',i2,',',i1,' =')
!!!!
!!!          do isph = 1,nsph
!!!            do i = 1,nbasis
!!!              write(*,"(4x,300(e12.5,2x))") ( dA((isph-1)*nbasis+i,j,ksph,icomp), j=1,nbasis*nsph )
!!!            enddo
!!!          enddo
!!!!          
!!!!         numerical derivatives
!!!          write(*,1006) ksph,icomp
!!! 1006     format(' ( A(r+r_',i2,',',i1,') - A(r) ) / eps =')
!!!! 
!!!          do isph = 1,nsph
!!!            do i = 1,nbasis
!!!              write(*,"(4x,300(e12.5,2x))") &
!!!              ( ( A_plus((isph-1)*nbasis+i,j,ksph,icomp)- &
!!!                       A((isph-1)*nbasis+i,j)             ) / eeps , j=1,nbasis*nsph )
!!!            enddo
!!!          enddo
!!!          write(*,*)''
!!!! 
!!!        enddo
!!!      enddo
!
!     restore DS
      call memfree
      call ddinit( nsph_save, x_save, y_save, z_save, r_save )
!
!     deactivate quiet flag      
      iquiet = .false.
!
!
endsubroutine ADJcheck

!---------------------------------------------------------------------------------------
subroutine check_forcesPCM( Psi0, sigma0, charge, f )
!
      use ddcosmo , only : nbasis, nsph, iquiet, csph, rsph, memfree, ddinit, &
                           eps, ncav, ccav, ngrid, zero, sprod, wghpot, one, &
                           lmax
!                           
      implicit none
      real*8, dimension(nbasis,nsph), intent(in) :: Psi0
      real*8, dimension(nbasis,nsph), intent(in) :: sigma0
      real*8, dimension(       nsph), intent(in) :: charge
      real*8, dimension(     nsph,3), intent(in) :: f
!
      integer,parameter :: niter = 4
!
      real*8 :: phi(ncav), psi(nbasis,nsph), g(ngrid,nsph), sigma(nbasis,nsph)
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), phi_eps(nbasis,nsph)
      real*8 :: x(nsph), y(nsph), z(nsph), rwork(niter,nsph*3), rrate(niter,nsph*3), &
                hwork(niter,nsph*3),ework(niter,nsph*3),fwork(niter,nsph*3)
      real*8 :: E0, E_plus, err, eeps, h
      integer :: iter, icomp, ksph, nsph_save, j, ncav_save, nbasis_save
!
      character(len=10) :: x1,x2
      character(len=30) :: fname
      integer :: fp
!
!---------------------------------------------------------------------------------------
!
!     activate quiet flag
      iquiet = .true.
!
!     compute E0
      E0 = 0.5d0 * (eps-1.d0)/eps * sprod( nbasis*nsph, sigma0, Psi0 )
!
!     save initial DS
      nsph_save = nsph
      ncav_save = ncav
      nbasis_save = nbasis
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
!     set initial increment
      eeps=0.01d0
!      
!     initialize
      rwork = zero ; rrate = zero ; ework = zero ; fwork = zero
!
!     loop over increments
      do iter = 1,niter
!        
!       loop over d / dr_k
        do ksph = 1,nsph_save
!
!         loop over components of d / dr_k
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph)*(1.d0+eeps) ; h = eeps*x_save(ksph)
            case(2) ; y(ksph) = y_save(ksph)*(1.d0+eeps) ; h = eeps*y_save(ksph)
            case(3) ; z(ksph) = z_save(ksph)*(1.d0+eeps) ; h = eeps*z_save(ksph)
            endselect
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!            
!           potential Phi and Psi vector
            call mkrhs( nsph_save, charge, x, y, z, ncav_save, ccav, phi, nbasis_save, psi )
!
!           solve PCM equations       
            g(:,:)=zero ; sigma(:,:)=zero
            call wghpot( phi, g )
            call iefpcm( g, psi, sigma, phi_eps )
!
!           compute energy
            E_plus = 0.5d0 * (eps-1.d0)/eps * sprod( nbasis*nsph, sigma, psi )
!
!           account for a potentially null shift
            if ( abs(h) .gt. 1.E-12 ) then
!
!             account for a potentially null component of the force
              if ( abs( f(ksph,icomp) ).gt.1.E-12 ) then
!                        
!               compute relative error
!!!                err = abs( (E_plus - E0) / eeps + f(ksph,icomp) ) / abs( f(ksph,icomp) )
                err = abs( (E_plus - E0) / h + f(ksph,icomp) ) / abs( f(ksph,icomp) )
!
!               store
                rwork(iter,(ksph-1)*3+icomp) = err
                hwork(iter,(ksph-1)*3+icomp) = h
                ework(iter,(ksph-1)*3+icomp) = ( E_plus - E0 ) / h
                fwork(iter,(ksph-1)*3+icomp) = f(ksph,icomp)

                
!!!                if ( abs( ( E_plus - E0 ) / h ) .gt. 1.E+03) then
!!!                  write(*,*) 'E_plus = ',E_plus
!!!                  write(*,*) 'E0     = ',E0
!!!                  stop
!!!                endif

!               compute rate
                if ( iter.gt.1 ) then 
                  rrate(iter,(ksph-1)*3+icomp) =  log( rwork(iter-1,(ksph-1)*3+icomp) / &
                                                       rwork(iter  ,(ksph-1)*3+icomp)   ) / &
                                                  log( hwork(iter  ,(ksph-1)*3+icomp) / &
                                                       hwork(iter-1,(ksph-1)*3+icomp)   )
                endif
              endif
            endif
!
          enddo
        enddo
!
        eeps = eeps / 2.d0
!
      enddo
!      
!     print numerical derivative of energy
      write(*,*)'Numerical derivative of energy : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( ework(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!      
!     print analytical force
      write(*,*)'Analytical force : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' force  _'i2','i1' : ',300(e12.5,2x))") j,icomp, ( fwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!     print relative error
      write(*,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!      
!     print rate of convergence
      write(*,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!     print to file
!     -------------
!
      write(x1,'(I2.2)') lmax
      write(x2,'(I4.4)') ngrid
      fname = 'lmax' // trim(x1) // '_ngrid' // trim(x2)
      fp    = 17
!
      open( unit=fp, file=fname, form='formatted', access='sequential', status='unknown')
!      
      write(fp,*)'lmax,ngrid = ',lmax,ngrid
      write(fp,*) ''
!
!     print numerical derivative of energy
      write(fp,*)'Numerical derivative of energy : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( ework(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!      
!     print analytical force
      write(fp,*)'Analytical force : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' force  _'i2','i1' : ',300(e12.5,2x))") j,icomp, ( fwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!
!     print relative error
      write(fp,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!      
!     print rate of convergence
      write(fp,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
!
      close(fp)
!      
!     restore DS
      call memfree
      call ddinit( nsph_save, x_save, y_save, z_save, r_save )
!
!     deactivate quiet flag      
      iquiet = .false.
!
!
endsubroutine check_forcesPCM
