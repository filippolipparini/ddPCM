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
! Since A_eps' is independent of eps, let just set A' = A_oo' = A_eps' , so that :
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
!     STEP 1 : solve adjoint problem (A_eps L)^T s = Psi
!     --------------------------------------------------
!   
!     solve L^T y = Psi     
      call itsolv2( star, .true., Psi, Psi, y, rvoid )
!
!     solve A_eps^T s = y
!      s = y
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
! compute f = f + < s , A' x >
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
      real*8 :: vvij,t_ijn,f3
      integer :: n,icomp,jcomp,jsph,ksph
!      
!------------------------------------------------------------------------------      
!
!     case j \ne i
!     ============
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
!             compute derivatives of t_ijn
              dt_ijn(1:3) = rsph(isph) * vij(1:3) / vvij**3
!
!             compute derivatives of s_ijn
              do icomp = 1,3
                do jcomp = 1,3
!
                  ds_ijn(icomp,jcomp) = - vij(icomp)*vij(jcomp) / vvij**3
!
                  if ( icomp.eq. jcomp ) then
!
                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) + one / vvij
!
                  endif
!                  
                enddo
              enddo
!             
!             compute Y_ll^mm(s_ijn)
              call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!             compute f2(j,n)
!             ---------------
              call compute_f2_jn( t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,jsph), f2(1:3) )
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
!     case j \eq i [1]
!     ================
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
!           compute s_ikn, t_ikn
            vij    = csph(:,isph) + rsph(isph)*grid(:,n) - csph(:,ksph)
            vvij   = sqrt(dot_product(vij,vij))
            t_ijn  = rsph(ksph)/vvij 
            s_ijn  =        vij/vvij
!
!           compute derivatives of t_ikn
            dt_ijn(1:3) = - rsph(ksph) * vij(1:3) / vvij**3
!
!           compute derivatives of s_ikn
            do icomp = 1,3
              do jcomp = 1,3
!
                ds_ijn(icomp,jcomp) = - vij(icomp)*vij(jcomp) / vvij**3
!
                if ( icomp.eq. jcomp ) then
!
                  ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) + one / vvij
!
                endif
!                
              enddo
            enddo
!           
!           compute Y_ll^mm(s_ikn)
            call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!           compute f2(j,n)
!           ---------------
            call compute_f2_kn( ksph, n, t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,ksph), f2(1:3) )
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
!     case j \eq i [2]
!     ================
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize f4(n)
        f4(1:3) = zero
!
!       compute f3(n) [ store in f2 for convenience ]
!       ---------------------------------------------
        f2(1:3) = zi(1:3,n,isph) * dot_product( basis(:,n), s(:,isph) )
!
!       compute f2(n) [ store in f3 for convenience ]
!       ---------------------------------------------
        call compute_f2_n( basloc, x(:,isph), f3 )
!                
!       accumulate for f4(n)
!       --------------------
        f4(1:3) = f4(1:3) + f2(1:3) * f3
!
!       accumulate for f
!       ----------------
        f(1:3) = f(1:3) + w(n) * f4(1:3)
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
!              = t^l [ (l+1) \grad t Y(s) + t ( \grad s )^T \grad Y(s) ]
!
!              = tt [ (l+1) s1(1:3,lm) + s2(1:3,lm) ]
!
! and :
!
!                  4 pi l
!   f2(1:3) = sum  ------  sum x(lm) * f(1:3,lm)
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
          s2(icomp,1:nbasis) = s2(icomp,1:nbasis) + dbsloc(jcomp,1:nbasis)*ds(icomp,jcomp)
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
        fac = four*pi*(fl+one)/(two*fl+one)
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
subroutine compute_f2_kn( ksph, n, t, dt, ds, basloc, dbsloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, four, pi, ui, zi
!
      implicit none
      integer,                     intent(in)  :: ksph
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
!              = ( \grad U ) t^l+1 * Y(s) + U * t^l [ (l+1) \grad t Y(s) + t ( \grad s )^T \grad Y(s) ]
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
        s1(icomp,1:nbasis) = dt(icomp)*basloc(1:nbasis)
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
          s2(icomp,1:nbasis) = s2(icomp,1:nbasis) + dbsloc(jcomp,1:nbasis)*ds(icomp,jcomp)
!
        enddo
!
        s2(icomp,1:nbasis) = ui(n,ksph)*s2(icomp,1:nbasis) + zi(icomp,n,ksph)*basloc(1:nbasis)
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
        fac = four*pi*(fl+one)/(two*fl+one)
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
!       build factor : 4pi*l / (2l+1)
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
        iflag=1
      endif

!     nothing for now ...   
      dphi(:,:,:,:) = 0.d0
!
!
endsubroutine compute_dphi
