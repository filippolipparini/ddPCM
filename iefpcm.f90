!-------------------------------------------------------------------------------
! Purpose : solve PCM equation R_eps S \sigma = -R_oo \Phi by solving the 
!           system :
!
!             { R_eps \Phi_eps = R_oo \Phi
!             {
!             {       S \sigma = -\Phi_eps    ( <== COSMO )
!
!           After numerical discretization, we obtain linear systems 
!
!             { A_eps W = A_oo g
!             {
!             {     L X = -W
!
! Arguments :
!
!   -sigma_g 
!-------------------------------------------------------------------------------
subroutine iefpcm( phi, psi, sigma_g )
!
      use  ddcosmo
!      
      implicit none
!
      real*8, dimension( ngrid,nsph), intent(in)    :: phi
      real*8, dimension(nbasis,nsph), intent(in)    :: psi
      real*8, dimension(nbasis,nsph), intent(inout) :: sigma_g
!      
!     P. Gatto, Nov 2016      
!     real*8, dimension(ngrid, nsph), intent(inout) :: sigma_g
!
!     local arrays:
      real*8, allocatable :: philm(:,:), wlm(:,:), glm(:,:), vold(:,:)
!
!     scratch arrays:
      real*8, allocatable :: basloc(:), vplm(:), vcos(:), vsin(:), xlm(:), x(:)
!
!     diis arrays:
      real*8, allocatable :: xdiis(:,:,:), ediis(:,:,:), bmat(:)
!
!     preconditioner:
      real*8, allocatable :: prec(:,:,:), precm1(:,:,:)
!
      integer :: it, isph, nmat, lenb, istatus
      real*8  :: ene, vrms, vmax, tol
      real*8, allocatable :: err(:), ddiag(:)
      logical :: dodiis
!
      real*8, allocatable :: phi_eps(:,:), dphi(:,:,:,:), f(:,:)
!
      integer, parameter :: nitmax=300
      real*8,  parameter :: ten=10.0d0, tredis=1.d-2
!
!-------------------------------------------------------------------------------
!
!     set solver tolerance
      tol = ten**(-iconv)
!      
!     set up DIIS iterative solver
      dodiis = .false.
      nmat   = 1
      lenb   = ndiis + 1
!
!     allocate memory
      allocate( err(   nbasis)          , &
                xlm(   nbasis)          , &
                x(     ngrid)           , &
                vold(  nbasis,nsph)     , &
                philm( nbasis,nsph)     , &
                wlm(   nbasis,nsph)     , &
                glm(   nbasis,nsph)     , &
                basloc(nbasis)          , &  ! spherical harmonics (SH)
                vplm(  nbasis)          , &  ! service array for SH
                vcos(  lmax+1)          , &  ! service array for SH
                vsin(  lmax+1)          , &  ! service array for SH
                xdiis(nbasis,nsph,ndiis), &
                ediis(nbasis,nsph,ndiis), &
                bmat(lenb*lenb)         , &
                prec(nbasis,nbasis,nsph), precm1(nbasis,nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : allocation failed !'
        stop
      endif
!
!
!===================================================================================
! PCM                                                                              |
!===================================================================================
!
!     STEP 1 : compute rhs phi = A_oo g
!     ---------------------------------
!
!     loop over atoms
      do isph = 1,nsph
!      
!       compute SH expansion glm of phi
        call intrhs( isph, phi(:,isph), glm(:,isph) )
!        
      end do
!
!     P. Gatto, Dec 2016 : why is it not initialized with philm ???
!
!     initial residual : R^0  = g - A_eps * x^0 = g
      vold(:,:) = glm(:,:)
!      
      !$omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!      
!     loop over atoms
      do isph = 1,nsph
!     
!       phi = A_oo g
        call mkrvec( isph, zero, glm, philm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!        
      end do
!
      !$omp parallel do default(shared) private(isph)
!      
!      
!     STEP 2 : solve A_eps W = phi
!     -----------------------------
!    
!     loop over atoms
      do isph = 1,nsph
!
!       compute inverse of diagonal block A_eps,ii
        call mkprec( isph, .true., prec(:,:,isph), precm1(:,:,isph) )
!        
      end do
!
!
!                                   n
!     STEP 2.1 : Jacobi method for W
!     -------------------------------
!
!     main loop
      write(iout,1000)
      write(iout,1010)
 1000 format('   first loop: computing V(eps)')
 1010 format('   it        error        err-00')
!      
!     Jacobi iteration
      do it = 1,nitmax
!      
!       initialize residual to zero
        wlm(:,:) = zero
!
      !$omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!
!                                    n
!       STEP 2.2 : compute residual R
!       ------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
          call mkrvec( isph, eps, vold, wlm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!          
        end do
!
!       R^n = Phi - A_eps * R^n-1 
        wlm(:,:) = philm(:,:) - wlm(:,:)
!
!                                n          
!       STEP 2.3 : solve for  W_i  = A_eps,ii^-1 ( ... )
!       ------------------------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
!         apply inverse
          call DGEMV( 'N', nbasis, nbasis, one, precm1(:,:,isph), nbasis, wlm(:,isph), 1, zero, xlm, 1 )
!
!         update : W_i^n = Phi_i - A_eps,ii W^n-1 
          wlm(:,isph) = vold(:,isph) + xlm(:)
!          
        end do
!
!               n-1    n
!       update W    = W
        vold = vold - wlm
!
!
!       STEP 3.3 : check for convergence
!       --------------------------------
!
        err(:) = zero
!
!       loop over atoms
        do isph = 1,nsph
!        
!         accumulate
          err(:) = err(:) + vold(:,isph)**2
!          
        end do
!
        err = sqrt(err/dble(nsph))
!
!       compute rms- and max-norm of v_old
        call rmsvec( nbasis*nsph, vold, vrms, vmax )
!
        if ( vrms.le.tredis .and. ndiis.gt.0 )  dodiis = .true.
!
!       diis extrapolation
        if ( dodiis ) then
!                
          xdiis(:,:,nmat) = wlm
          ediis(:,:,nmat) = vold
!          
          call diis( nbasis*nsph, nmat, xdiis, ediis, bmat, wlm )
!          
        end if
!        
        write(iout,1020) it, vrms, err(1)
 1020   format(1x,i4,f14.8,121d12.4)
!
!       convergence has been achieved
        if ( vrms.lt.tol )  goto 999
!        
!       update
        vold = wlm
!        
      end do
!
!     convergence has not been reached
      stop ' convergence failure!!!'
!
  999 continue
!  
      write(iout,2000)
 2000 format('   first loop has converged.',/,'   second loop: solving ddCOSMO equations for V(eps)')
!
!     compute charge distribution and energy
!
!===================================================================================
! ddCOSMO                                                                          |
!===================================================================================
!
!     initialize
      sigma_g = zero
!
!     solve  L sigma = W , compute energy
      call itsolv2( .false., .true., wlm, psi, sigma_g, ene )
!
!   
!===================================================================================
! FORCES                                                                           |
!===================================================================================
!
!     allocate 
      allocate( phi_eps(nbasis,nsph), dphi(ngrid,nsph,nsph,3), f(3,nsph) , stat=istatus )
!      
!     check
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : [1] failed allocation !'
        stop
      endif 
!
!     redirect
      phi_eps = -wlm
!      
!     provide derivatives of Phi
      call compute_dphi( dphi )
!      
!     compute forces
      call compute_forces( phi, dphi, psi, sigma_g, phi_eps, f )
!      
!     deallocate
      deallocate( phi_eps, dphi, f , stat=istatus )
!      
!     check
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : [1] failed deallocation !'
        stop
      endif 
!
!
!     free the memory
      deallocate( err, xlm, x, vold, philm, wlm, glm, basloc, vplm, vcos, vsin, &
                  xdiis, ediis, bmat, prec, precm1 , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : deallocation failed !'
        stop
      endif
!      
      return
!
!
endsubroutine iefpcm
!----------------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! Purpose : solve PCM equation A_eps W = Phi
!-------------------------------------------------------------------------------
subroutine ADJpcm( philm, wlm )
!
      use  ddcosmo
!      
      implicit none
!
      real*8, dimension(nbasis,nsph), intent(in ) :: philm
      real*8, dimension(nbasis,nsph), intent(out) :: wlm
!      
!     local arrays:
      real*8, allocatable :: vold(:,:)
!
!     scratch arrays:
      real*8, allocatable :: basloc(:), vplm(:), vcos(:), vsin(:), xlm(:), x(:)
!
!     diis arrays:
      real*8, allocatable :: xdiis(:,:,:), ediis(:,:,:), bmat(:)
!
!     preconditioner:
      real*8, allocatable :: prec(:,:,:), precm1(:,:,:)
!
      integer :: it, isph, nmat, lenb, istatus
      real*8  :: ene, vrms, vmax, tol
      real*8, allocatable :: err(:), ddiag(:)
      logical :: dodiis
!
      integer, parameter :: nitmax=300
      real*8,  parameter :: ten=10.0d0, tredis=1.d-2
!
!-------------------------------------------------------------------------------
!
!     set solver tolerance
      tol = ten**(-iconv)
!      
!     set up DIIS iterative solver
      dodiis = .false.
      nmat   = 1
      lenb   = ndiis + 1
!
!     allocate memory
      allocate( err(   nbasis)          , &
                xlm(   nbasis)          , &
                x(     ngrid)           , &
                vold(  nbasis,nsph)     , &
                basloc(nbasis)          , &  ! spherical harmonics (SH)
                vplm(  nbasis)          , &  ! service array for SH
                vcos(  lmax+1)          , &  ! service array for SH
                vsin(  lmax+1)          , &  ! service array for SH
                xdiis(nbasis,nsph,ndiis), &
                ediis(nbasis,nsph,ndiis), &
                bmat(lenb*lenb)         , &
                prec(nbasis,nbasis,nsph), precm1(nbasis,nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'ADJpcm : allocation failed !'
        stop
      endif
!
!     initial residual : R^0  = phi - A_eps * x^0 = phi
      vold(:,:) = philm(:,:)
!      
!     STEP 2 : solve A_eps W = phi
!     -----------------------------
!    
!     loop over atoms
      do isph = 1,nsph
!
!       compute inverse of diagonal block A_eps^T,ii
        call ADJprec( isph, .true., prec(:,:,isph), precm1(:,:,isph) )
!        
      enddo
!
!
!                                   n
!     STEP 2.1 : Jacobi method for W
!     -------------------------------
!
!     main loop
      write(iout,1000)
      write(iout,1010)
 1000 format('   first loop: computing V(eps)')
 1010 format('   it        error        err-00')
!      
!     Jacobi iteration
      do it = 1,nitmax
!      
!       initialize residual to zero
        wlm(:,:) = zero
!
      !$omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!
!                                    n
!       STEP 2.2 : compute residual R
!       ------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
          call ADJvec( isph, eps, vold, wlm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!          
        end do
!
!       R^n = Phi - A_eps * R^n-1 
        wlm(:,:) = philm(:,:) - wlm(:,:)
!
!                                n          
!       STEP 2.3 : solve for  W_i  = A_eps,ii^-1 ( ... )
!       ------------------------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
!         apply inverse
          call DGEMV( 'N', nbasis, nbasis, one, precm1(:,:,isph), nbasis, wlm(:,isph), 1, zero, xlm, 1 )
!
!         update : W_i^n = Phi_i - A_eps,ii W^n-1 
          wlm(:,isph) = vold(:,isph) + xlm(:)
!          
        end do
!
!               n-1    n
!       update W    = W
        vold = vold - wlm
!
!
!       STEP 3.3 : check for convergence
!       --------------------------------
!
        err(:) = zero
!
!       loop over atoms
        do isph = 1,nsph
!        
!         accumulate
          err(:) = err(:) + vold(:,isph)**2
!          
        end do
!
        err = sqrt(err/dble(nsph))
!
!       compute rms- and max-norm of v_old
        call rmsvec( nbasis*nsph, vold, vrms, vmax )
!
        if ( vrms.le.tredis .and. ndiis.gt.0 )  dodiis = .true.
!
!       diis extrapolation
        if ( dodiis ) then
!                
          xdiis(:,:,nmat) = wlm
          ediis(:,:,nmat) = vold
!          
          call diis( nbasis*nsph, nmat, xdiis, ediis, bmat, wlm )
!          
        end if
!        
        write(iout,1020) it, vrms, err(1)
 1020   format(1x,i4,f14.8,121d12.4)
!
!       convergence has been achieved
        if ( vrms.lt.tol )  goto 999
!        
!       update
        vold = wlm
!        
      end do
!
!     convergence has not been reached
      stop ' convergence failure!!!'
!
  999 continue
!  
      write(iout,2000)
 2000 format('   first loop has converged.',/,'   second loop: solving ddCOSMO equations for V(eps)')
!
!
!
!     free the memory
      deallocate( err, xlm, x, vold, basloc, vplm, vcos, vsin, &
                  xdiis, ediis, bmat, prec, precm1 , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'ADJpcm : deallocation failed !'
        stop
      endif
!
!
endsubroutine ADJpcm
!----------------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------------
! Purpose :
!
! Compute on-the-fly action of operator :
!
!   dvlm_i = A_eps,i * vlm 
!
!          = sum A_eps,ij * vlm_j
!             j
!                                                            ~
!          = 2pi (eps+1) / (eps-1) vlm_i - D_i vlm_i -   sum D_ij vlm_j
!                                                      j \ne i
!
! where :
!
!                        2pi
!   D_i vlm_i = -  sum  -----  sum w_n  Y_l^m(s_n)  U_i^n  Y_l'^m'(s_n)   vl'm'_i
!                 l',m' 2l'+1   n      
!
!                                                   -2pi
!             = sum  w_n  Y_l^m(s_n)  U_i^n   sum  -----  Y_l'^m'(s_n)  vl'm'_i
!                n                           l',m' 2l'+1
!
!                                     |---------------- x(i,n) ---------------|
!
!
!        ~                            4pi l'                                     l'+1
!   sum  D_ij vlm_j = -  sum     sum  ----- sum w_n  Y_l^m(s_n)  U_j^n  ( t_ijn )      Y_l'^m' (s_ijn)  vl'm'_j
! j \ne i              j \ne i  l',m' 2l'+1  n
!   
!                                                                -4pi l'           l'+1
!                   = sum  w_n  Y_l^m(s_n)  U_j^n   sum     sum  -------  ( t_ijn )      Y_l'^m'(s_ijn)  vl'm'_j
!                      n                          j \ne i  l',m'  2l'+1
! 
!                                           |---------------------------- x(i,n) ------------------------------|
!
!
! Remark : when eps_s=0, the eps=oo case is triggered, i.e., 2pi Id - D_i - ...
!----------------------------------------------------------------------------------------
!
subroutine mkrvec( isph, eps_s, vlm, dvlm, xlm, x, basloc, vplm, vcos, vsin )
!
      use  ddcosmo
!     use  dd_utils , only : nbasis, nsph, ngrid, lmax, csph, rsph, grid, basis, ui, facl
!      
      implicit none
      integer,                         intent(in   ) :: isph
      real*8,                          intent(in   ) :: eps_s
      real*8,  dimension(nbasis,nsph), intent(in   ) :: vlm
      real*8,  dimension(nbasis),      intent(inout) :: dvlm
      real*8,  dimension(ngrid),       intent(inout) :: x
      real*8,  dimension(nbasis),      intent(inout) :: xlm, basloc, vplm
      real*8,  dimension(lmax+1),      intent(inout) :: vcos, vsin
!
      integer :: its, jsph
      real*8  :: vij(3), sij(3), fep
      real*8  :: vvij, tij, stslm, stslm2, stslm3
!
!----------------------------------------------------------------------------------------
!
!     compute f( \eps )
      fep = two*pi*(eps_s+one)/(eps_s-one)
      if ( eps_s.eq.zero )  fep = two*pi
!
!     initialize
      x(:) = zero
!
!     loop over grid points of i-sphere 
      do its = 1,ngrid
! 
!       non-null contribution from integration point 
        if ( ui(its,isph).gt.zero ) then
!
!         loop over spheres
          do jsph = 1, nsph
!
!           action of D~_ij
!           ---------------
!
            if ( jsph.ne.isph ) then
!
!             compute t_ijn, s_ijn
              vij  = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
              vvij = sqrt(dot_product(vij,vij))
              tij  = vvij/rsph(jsph) 
              sij  = vij/vvij
!              
!             compute Y_l'^m'(s_ijn)
              call ylmbas( sij, basloc, vplm, vcos, vsin )
!              
!             accumulate for x
              if ( tij.lt.one ) then
!                      
                x(its) = x(its) + dtslm(  tij, vlm(:,jsph), basloc )
!                      
              else
!                      
                x(its) = x(its) + dtslm2( tij, vlm(:,jsph), basloc )
!                      
              end if
!                      
!                      
!           action of D_ii
!           --------------
!                      
            else
!
              xlm(:) = -pt5 * vlm(:,isph) / facl(:)
!
!             accumulate for x
              x(its) = x(its) + dot_product( basis(:,its), xlm(:) )
!              
            end if
          end do
        end if
      end do
!
!     compute x
      x(:) = x(:) * ui(:,isph)
!
!     integrate against SH
      call intrhs( isph, x, dvlm )
!
!
!     add action of 2pi * f(eps) * Id
!     -------------------------------
!
      dvlm(:) = fep * vlm(:,isph) - dvlm(:)
!
      return
!
!
endsubroutine mkrvec
!-------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------------
! Purpose :
!
! Compute on-the-fly action of operator :
!
!                 T
!   dvlm_i = A_eps ,i * vlm 
!
!                     T
!          = sum A_eps ,ij * vlm_j
!             j
!                                                            ~
!          = 2pi (eps+1) / (eps-1) vlm_i - D_i vlm_i -   sum D_ij vlm_j
!                                                      j \ne i
!
! where :
!
!                 2 pi
!   D_i vlm_i = - ----  sum  w_n  U_i^n  Y_l^m(s_n)   sum   Y_l'^m'(s_n)  vl'm'_i
!                 2l+1   n                           l',m'
!
! and :
!
!        ~              4pi l                                                   l+1
!   sum  D_ij vlm_j = - -----  sum  w_n    sum    Y_l^m(s_ijn)  U_j^n  ( t_ijn )      sum   Y_l'^m'(s_n)  vl'm'_j
! j \ne i                2l+1   n        j \ne i                                     l',m'
!   
!
! Remark : when eps_s=0, the eps=oo case is triggered, i.e., 2pi Id - D_i - ...
!----------------------------------------------------------------------------------------
!
subroutine ADJvec( isph, eps_s, vlm, dvlm, xlm, x, basloc, vplm, vcos, vsin )
!
      use  ddcosmo , only : nbasis, nsph, ngrid, lmax, csph, rsph, grid, basis, ui, facl, &
                            facll, two, pi, one, zero, pt5, w, ylmbas
!      
      implicit none
      integer,                         intent(in   ) :: isph
      real*8,                          intent(in   ) :: eps_s
      real*8,  dimension(nbasis,nsph), intent(in   ) :: vlm
      real*8,  dimension(nbasis),      intent(inout) :: dvlm
      real*8,  dimension(ngrid),       intent(inout) :: x
      real*8,  dimension(nbasis),      intent(inout) :: xlm, basloc, vplm
      real*8,  dimension(lmax+1),      intent(inout) :: vcos, vsin
!
      integer :: n, jsph, l, m, ind
      real*8 :: vij(3), s_ijn(3), fep
      real*8 :: vvij, t_ijn, tt, ss
      real*8 :: dij_vlm(nbasis), dii_vlm(nbasis), di_vlm(nbasis), dj_vlm(nbasis), &
                f1(nbasis)
!
!----------------------------------------------------------------------------------------
!
!     compute f( \eps )
      fep = two*pi*(eps_s+one)/(eps_s-one)
      if ( eps_s.eq.zero )  fep = two*pi
!
!     initialize [ this is just a DUMMY variable ]
      x(:) = zero
!
!     loop over grid points of i-sphere 
      do n = 1,ngrid
!
!       initialize
        dij_vlm(:)=zero ; dii_vlm(:)=zero
!
!       loop over spheres
        do jsph = 1,nsph
!          
!         non-null contribution from integration point 
          if ( ui(n,jsph).gt.zero ) then
!
!
!           action of D~_ij
!           ---------------
!
            if ( jsph.ne.isph ) then
!
!             compute t_ijn, s_ijn
              vij   = csph(:,jsph) + rsph(jsph)*grid(:,n) - csph(:,isph)
              vvij  = sqrt( dot_product( vij,vij ) )
              t_ijn = rsph(isph)/vvij 
              s_ijn =        vij/vvij
!              
!             contract over ll, mm
              ss = dot_product( basis(:,n), vlm(:,jsph) )
!              
!             compute Y_l^m(s_ijn)
              call ylmbas( s_ijn, basloc, vplm, vcos, vsin )
!              
!             compute f1(lm) = t_ijn^l+1 Y_l^m(s_ijn)
              if ( t_ijn.lt.one ) then

!               initialize t_ijn^l+1 factor
                tt = one
!
                do l = 0,lmax
!                
!                 update factor
                  tt = tt*t_ijn
!                  
                  ind = l*l + l + 1
!                  
                  do m = -l,l
!                   
                    f1(ind+m) = tt * basloc(ind+m)
!                    
                  enddo
                enddo
              endif
!
!             compute f1(lm) = U_j^n t_ijn^l+1 Y_l^m(s_ijn) * ss
              f1(:) = ui(n,jsph) * f1(:) * ss
!
!             accumulate over j
              dij_vlm(:) = dij_vlm(:) + f1(:)
!                      
!                      
!           action of D_ii
!           --------------
!                      
            else
!
!             contratct over ll, mm   
              ss = dot_product( basis(:,n), vlm(:,isph) )
!              
!             compute U_i^n Y_l^m(s_n) * ss              
              dii_vlm(:) = ui(n,isph) * basis(:,n) * ss
!              
            endif
!
          endif
        enddo
!        
!       accumulate over n
        dj_vlm(:) = dj_vlm(:) + w(n) * dij_vlm(:)
        di_vlm(:) = di_vlm(:) + w(n) * dii_vlm(:)
!
      enddo
!
!
!     add action of 2pi * f(eps) * Id
!     -------------------------------
!
      dvlm(:) = fep*vlm(:,isph) - pt5/facl(:)*di_vlm(:) - facll(:)*dj_vlm(:)
!
      return
!
!
endsubroutine ADJvec
!-------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! Purpose : compute A_eps,ii^-1.
!
! Remark : this is a bit more than a preconditioner... in fact, it is the
!          exact (numerical) inverse. 
!
!-------------------------------------------------------------------------------
!
subroutine mkprec( isph, doinv, prec, precm1 )
!
      use  ddcosmo
!     use  dd_utils , only : nbasis, ngrid, basis, lmax, ui, w, eps
!      
      implicit none
      integer,                          intent(in   ) :: isph
      logical,                          intent(in   ) :: doinv
      real*8, dimension(nbasis,nbasis), intent(inout) :: prec
      real*8, dimension(nbasis,nbasis), intent(inout) :: precm1
!
      integer :: l, m, ind, l1, m1, ind1, info, its, istatus
      real*8  :: f, f1
!
      integer, allocatable :: ipiv(:)
      real*8,  allocatable :: work(:)
!
!-------------------------------------------------------------------------------
!
      allocate( ipiv(nbasis) , work(nbasis*nbasis) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'mkprec : allocation failed !'
        stop
      endif
!
!     initialize
      prec(:,:) = zero
!      
!      
!     STEP 1
!     ------
!
!     loop over grid points
      do its = 1,ngrid
!      
!       1st loop over SH degree
        do l = 0,lmax
!
!         index associated to Y_l^0
          ind = l*l + l + 1
!
!         1st loop over SH order
          do m = -l,l
!          
            f = w(its) * basis(ind+m,its) * ui(its,isph)
!            
!           2nd loop over SH degree
            do l1 = 0,lmax
!
!             index associated to Y_l1^0
              ind1 = l1*l1 + l1 + 1
!
!             2nd loop over SH order
              do m1 = -l1, l1
!              
                f1 = two * pi / (two*dble(l1) + one) * basis(ind1+m1,its)
!                
!               accumulate
                prec(ind+m,ind1+m1) = prec(ind+m,ind1+m1) + f * f1
!                
              end do
            end do
          end do
        end do
      end do
!
!
!     STEP 2 : diagonal part
!     ----------------------
!
!     lopp over SH degree
      do l = 0,lmax
!
!       index associated to Y_l^0
        ind = l*l + l + 1
!
!       loop over SH order
        do m = -l,l
!        
          f = two * pi * (eps+one) / (eps-one)
!
!         accumulate
          prec(ind+m,ind+m) = prec(ind+m,ind+m) + f
!          
        end do
      end do
!      
!
!     STEP 3 : invert preconditioner
!     ------------------------------
!
      if ( doinv ) then
!              
        precm1(:,:) = prec(:,:)
!        
!       compute LU factorization
        call DGETRF( nbasis, nbasis, precm1, nbasis, ipiv, info )
!
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRF : info = ', info
        end if
!        
!       invert factorization
        call DGETRI( nbasis, precm1, nbasis, ipiv, work, nbasis*nbasis, info )
!        
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRI : info = ', info
        end if
!
      end if
!      
      deallocate ( work, ipiv, stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'mkprec : deallocation failed !'
        stop
      endif
!      
      return
!      
!      
endsubroutine mkprec
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
subroutine ADJprec( isph, doinv, prec, precm1 )
!
      use  ddcosmo
!     use  dd_utils , only : nbasis, ngrid, basis, lmax, ui, w, eps
!      
      implicit none
      integer,                          intent(in   ) :: isph
      logical,                          intent(in   ) :: doinv
      real*8, dimension(nbasis,nbasis), intent(inout) :: prec
      real*8, dimension(nbasis,nbasis), intent(inout) :: precm1
!
      integer :: l, m, ind, l1, m1, ind1, info, its, istatus
      real*8  :: f, f1
!
      integer, allocatable :: ipiv(:)
      real*8,  allocatable :: work(:)
!
!-------------------------------------------------------------------------------
!
      allocate( ipiv(nbasis) , work(nbasis*nbasis) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'ADJprec : allocation failed !'
        stop
      endif
!
!     initialize
      prec(:,:) = zero
!      
!      
!     STEP 1
!     ------
!
!     loop over grid points
      do its = 1,ngrid
!      
!       1st loop over SH degree
        do l = 0,lmax
!
!         index associated to Y_l^0
          ind = l*l + l + 1
!
!         1st loop over SH order
          do m = -l,l
!          
            f = two*pi/(two*dble(l) + one) * w(its) * basis(ind+m,its) * ui(its,isph)
!            
!           2nd loop over SH degree
            do l1 = 0,lmax
!
!             index associated to Y_l1^0
              ind1 = l1*l1 + l1 + 1
!
!             2nd loop over SH order
              do m1 = -l1, l1
!                
!               accumulate
                prec(ind1+m1,ind+m) = prec(ind1+m1,ind+m) + f * basis(ind1+m1,its)
!                
              end do
            end do
          end do
        end do
      end do
!
!
!     STEP 2 : diagonal part
!     ----------------------
!
!     loop over SH degree
      do l = 0,lmax
!
!       index associated to Y_l^0
        ind = l*l + l + 1
!
!       loop over SH order
        do m = -l,l
!        
          f = two * pi * (eps+one) / (eps-one)
!
!         accumulate
          prec(ind+m,ind+m) = -prec(ind+m,ind+m) + f
!          
        end do
      end do
!      
!
!     STEP 3 : invert preconditioner
!     ------------------------------
!
      if ( doinv ) then
!              
        precm1(:,:) = prec(:,:)
!        
!       compute LU factorization
        call DGETRF( nbasis, nbasis, precm1, nbasis, ipiv, info )
!
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRF : info = ', info
        end if
!        
!       invert factorization
        call DGETRI( nbasis, precm1, nbasis, ipiv, work, nbasis*nbasis, info )
!        
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRI : info = ', info
        end if
!
      end if
!      
      deallocate ( work, ipiv, stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'mkprec : deallocation failed !'
        stop
      endif
!      
      return
!      
!      
endsubroutine ADJprec
