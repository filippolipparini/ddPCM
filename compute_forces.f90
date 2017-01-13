!-------------------------------------------------------------------------------------
! Recall :
!
!   A_eps Phi_eps = A_oo Phi     <==  PCM
!        
!         L sigma = -Phi_eps     <==  COSMO
!      
! Computation of forces :      
!      
!   F = - < s , -A_oo' Phi - A_oo Phi' + A_eps' Phi_eps - A_eps L' sigma >
!
! where s is the solution of the adjoint problem 
!
!   (A_eps L)^T s = Psi
!
! Since A_eps' is independent of eps, let's just set A' = A_oo' = A_eps' , so that :
!
!   F = < s , A' ( Phi - Phi_eps ) >  + < s , A_oo Phi' > + < s , A_eps L' sigma >
!
!       < s , A' ( Phi - Phi_eps ) >  + < s , A_oo Phi' > + < A_eps^T s , L' sigma >
!
!-------------------------------------------------------------------------------------
!
subroutine compute_forces( Phi, dPhi, Psi, sigma, Phi_eps, f )
!
      use ddcosmo , only : zero, ngrid, nsph, nbasis, zero, lmax, intrhs, itsolv2, &
                           basis, fdoka, fdokb
!      
      implicit none
      real*8, dimension( ngrid,nsph),        intent(in)  :: Phi
      real*8, dimension( ngrid,nsph,nsph,3), intent(in)  :: dPhi
      real*8, dimension(nbasis,nsph),        intent(in)  :: Psi
      real*8, dimension(nbasis,nsph),        intent(in)  :: sigma
      real*8, dimension(nbasis,nsph),        intent(in)  :: Phi_eps
      real*8, dimension(nsph,3),             intent(out) :: f
!
      real*8, dimension(ngrid,nsph) :: xi
      real*8, dimension(nbasis,nsph) :: w_lm, s, y
      real*8, dimension(nbasis,nsph,nsph,3) :: g_lm
      real*8, dimension(nbasis,nsph,nsph,3) :: dphi_lm
      real*8, dimension(ngrid) :: x
      real*8, dimension(nbasis) :: xlm, basloc, vplm
      real*8, dimension(3,nbasis) :: dbsloc
      real*8, dimension(lmax+1) :: vcos, vsin
      real*8 :: rvoid
!
      integer :: isph,jsph,icomp,n
      logical, parameter :: star=.true.
!
!-------------------------------------------------------------------------------------
!
      call ADJcheck


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
!     STEP 2 : compute f = < s , A' ( Phi - Phi_eps ) >
!     -------------------------------------------------
!
!     initialize
      f(:,:) = zero
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
!
!     STEP 3 : compute f = f + < s , A_oo Phi' >
!     ------------------------------------------
!
!     compute SH expansion on i-sphere ...
      do isph = 1,nsph
!
!       ... of derivative of Phi with respect to center of j-sphere
        do jsph = 1,nsph
!       
          do icomp = 1,3
!       
            call intrhs( isph, dPhi(:,isph,jsph,icomp), g_lm(:,isph,jsph,icomp) )
!       
          enddo
        enddo
      enddo
!
!     apply A_oo
      do isph = 1,nsph
!
        do jsph = 1,nsph
!     
          do icomp = 1,3
!     
            call mkrvec( isph, zero, g_lm(:,:,jsph,icomp), dphi_lm(:,isph,jsph,icomp), xlm, x, basloc, vplm, vcos, vsin )
!
!           accumulate
            f(isph,icomp) = f(isph,icomp) + dot_product( s(:,isph) , dphi_lm(:,isph,jsph,icomp) )
!     
          enddo
        enddo
      enddo
!
!
!     STEP 4 : compute f = f - < y , L' sigma >
!     -----------------------------------------
!
!     compute xi
      do isph = 1,nsph
        do n = 1,ngrid
!        
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
                           ddinit,ui,zi
!
      implicit none
      real*8 :: f(nbasis,nsph),Af( nbasis)
      real*8 :: e(nbasis,nsph),ATe(nbasis)
      real*8 :: A(nbasis*nsph,nbasis*nsph),AT(nbasis*nsph,nbasis*nsph)
      real*8 :: dA1(nbasis*nsph,nbasis*nsph),dA2(nbasis*nsph,nbasis*nsph)
      real*8 :: dA3(nbasis*nsph,nbasis*nsph)
      real*8 :: dA(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: A_plus(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: xlm(nbasis),xx(ngrid),vplm(nbasis),vcos(lmax+1),vsin(lmax+1), &
                basloc(nbasis)
      integer :: isph,jsph,i,j,ibeg,iend,nsph_save,icomp,ksph,iter,n
      real*8 :: eps,s1,s2,eeps,err,rnorm
      integer, parameter :: niter = 6
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), s3(3)
      real*8 :: x(nsph), y(nsph), z(nsph), iwork(nsph*3,2), rwork(niter,nsph*3)
!
!--------------------------------------------------------------------------------
!
!     print A, A^T
!     ------------
      eps=two
!
      do isph = 1,nsph
!
        do jsph = 1,nsph
          do j = 1,nbasis
!
!           set basis vector
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
      write(*,*) 'A = '
      do isph = 1,nsph
        do i = 1,nbasis
          write(*,"(4x,300(e12.5,2x))") ( A((isph-1)*nbasis+i,j), j=1,nbasis*nsph )
        enddo
      enddo
      write(*,*)''
!
      write(*,*) '(A^T)^T = '
      do isph = 1,nsph
        do i = 1,nbasis
          write(*,"(4x,300(e12.5,2x))") ( AT(j,(isph-1)*nbasis+i), j=1,nbasis*nsph )
        enddo
      enddo
      write(*,*)''
!
!
!     check < A^T e , f > = < e , A f >
!     ---------------------------------
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
      e(:,:)=e(:,:)/sqrt(s1)
      f(:,:)=f(:,:)/sqrt(s2)
!
!     initialize
      Af( :)=zero
      ATe(:)=zero
      s1=zero
      s2=zero
!      
      do isph = 1,nsph
!
!       compute A_i f 
        call mkrvec( isph, eps, f, Af( :), xlm, xx, basloc, vplm, vcos, vsin )
!
!       compute A^T_i e
        call ADJvec( isph, eps, e, ATe(:), xlm, xx, basloc, vplm, vcos, vsin )
!
!       accumulate < e, A f >
        s1 = s1 + dot_product( e(:,isph), Af( :) )
!        
!       accumulate < A^T e, f >
        s2 = s2 + dot_product( f(:,isph), ATe(:) )
!        
      enddo
!      
      write(*,1000) abs(s1-s2) / abs(s1)
 1000 format(' | <e,Af> - <A^T f,e> | / |<e,Af>| = ', e12.5)
      write(*,*) ''
!
!!!      do isph = 1,nsph
!!!        do i = 1,ngrid
!!!          write(*,*)'n, isph, zi(;,n,isph) = ',i,isph,zi(:,i,isph)
!!!        enddo
!!!      enddo
!
!!!      do n=1,ngrid
!!!        do i=1,nsph
!!!          write(*,*)'n,i,ui,zi = ',n,i,ui(n,i),zi(1:3,n,i)
!!!        enddo
!!!      enddo
!
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
                e(:,:)=zero
                f(:,:)=zero
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
!     set increment initial
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
!           store error
            rwork(iter,(ksph-1)*3+icomp) = sqrt( err / rnorm )
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
!     printing
      do iter = 1,niter
!
        write(*,"(' iter = ',i1'; ',300(e12.5,2x))") iter, ( rwork(iter,j) , j=1,3*nsph)
!        
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
!
endsubroutine ADJcheck
