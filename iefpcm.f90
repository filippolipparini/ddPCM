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
!     STEP 1 : compute rhs R_oo Phi
!     -----------------------------
!
!     loop over atoms
      do isph = 1,nsph
!      
!       compute SH expansion of phi
        call intrhs( isph, phi(:,isph), glm(:,isph) )
!        
      end do
!
!     initial guess for solver
      vold(:,:) = glm(:,:)
!      
      !$omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!      
!     loop over atoms
      do isph = 1,nsph
!     
!       phi_l^m = R_oo g_l^m
        call mkrvec( isph, zero, glm, philm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!        
      end do
!
      !$omp parallel do default(shared) private(isph)
!      
!      
!     STEP 2 : build preconditioner
!     -----------------------------
!    
!     loop over atoms
      do isph = 1,nsph
!
        call mkprec( isph, .true., prec(:,:,isph), precm1(:,:,isph) )
!        
      end do
!
!
!     STEP 3 : solve R_eps Phi_eps = R_oo Phi 
!     ---------------------------------------
!
!     main loop
      write(iout,1000)
      write(iout,1010)
 1000 format('   first loop: computing V(eps)')
 1010 format('   it        error        err-00')
!      
!     loop until the max number of iterations
      do it = 1,nitmax
!      
!       initialize to zero
        wlm(:,:) = zero
!
      !$omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!
!       STEP 3.1 : compute action of R_eps
!       ----------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
!         w_l^m = R_eps v_old
          call mkrvec( isph, eps, vold, wlm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!          
        end do
!
!       compute b - A*x
        wlm(:,:) = philm(:,:) - wlm(:,:)
!
!       STEP 3.2 : apply preconditioner
!       -------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
!         x_l^m = precm1 * w_l^m
          call DGEMV( 'N', nbasis, nbasis, one, precm1(:,:,isph), nbasis, wlm(:,isph), 1, zero, xlm, 1 )
!
!         update w_l^m
          wlm(:,isph) = vold(:,isph) + xlm(:)
!          
        end do
!
!       update v_old
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
!     call to ddCOSMO
      sigma_g = zero
      call itsolv2( .false., .true., wlm, psi, sigma_g, ene )
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
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
! Purpose : compute action of operator
!
!             2pi (eps+1) / (eps-1) Id - sum D_ji
!                                         j
!
!           on the fly.
!
! Remark : when eps_s=0, the eps=oo case is triggered, i.e.,
!
!            2pi Id - D
!
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
!     loop over grid points
      do its = 1, ngrid
! 
!       positive contribution from integration point
        if ( ui(its,isph).gt.zero ) then
!
!         loop over spheres
          do jsph = 1, nsph
!
!
!           action of D_ji
!           --------------
!
            if ( jsph.ne.isph ) then
!
!             vv_n^jk
              vij  = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
!
!             v_n^jk
              vvij = sqrt(dot_product(vij,vij))
!
!             t_n^jk
              tij  = vvij/rsph(jsph) 
!
!             ss_n^jk
              sij  = vij/vvij
!              
!             compute spherical harmonics "basloc" at "sij"
              call ylmbas( sij, basloc, vplm, vcos, vsin )
!              
              if (tij.lt.one) then
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

              xlm(:) = -pt5 * vlm(:,isph) / facl(:)
              x(its) = x(its) + dot_product( basis(:,its), xlm(:) )
!              
            end if
          end do
        end if
      end do
!
      x(:) = x(:) * ui(:,isph)
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
