!-----------------------------------------------------------------------------------------
! Purpose : compute charge "sigma" and energy "ene" 
!-----------------------------------------------------------------------------------------
!
subroutine itsolv2( star, iefpcm, phi, psi, sigma, ene )
!
      use  ddcosmo
!     use  dd_utils , only : nbasis, nsph, eps, lmax, basis, eta, fact, csph, &
!                            fi, grid, iconv, iout, inl, nl, iprint, ngrid, w, &
!                            rsph, ndiis

!
      implicit none
      logical,                        intent(in   ) :: star
      logical,                        intent(in   ) :: iefpcm
      real*8, dimension(nbasis,nsph), intent(in   ) :: phi
      real*8, dimension(nbasis,nsph), intent(in   ) :: psi
      real*8, dimension(nbasis,nsph), intent(inout) :: sigma
      real*8,                         intent(inout) :: ene
!
!     local arrays
      real*8, allocatable :: pot(:), sigold(:,:), vlm(:)
!
!     scratch arrays
      real*8, allocatable :: basloc(:), vplm(:), vcos(:), vsin(:)
      real*8, allocatable :: delta(:), norm(:)
!
!     diis arrays
      real*8, allocatable :: xdiis(:,:,:), ediis(:,:,:), bmat(:), scr(:)
!
!     local variables
      integer :: it, isph, nmat, lenb, istatus
      real*8  :: tol, drms, dmax, sprod, fep
      logical :: dodiis, first
!
      integer, parameter :: nitmax=300
      real*8,  parameter :: ten=10.d0, tredis=1.0d-2
!
!-----------------------------------------------------------------------------------------
!
!     compute f(eps)
      fep = one
      if ( .not.iefpcm )  fep = (eps-one)/eps
!
!     set tolerance for solver
      tol = ten**(-iconv)
!
!     set up DIIS
      dodiis = .false.
      nmat   = 1
      lenb   = ndiis + 1
!
!     allocate workspaces
      allocate( pot(ngrid), vlm(nbasis), sigold(nbasis,nsph), delta(nbasis), &
                norm(nsph), vplm(nbasis), basloc(nbasis), vcos(lmax+1), &
                vsin(lmax+1) , xdiis(nbasis,nsph,ndiis), ediis(nbasis,nsph,ndiis), &
                bmat(lenb*lenb), scr(2*lenb*lenb) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolve2 : allocation failed !'      
        stop
      endif
!
!     initialize sigma^Jac,n-1
      sigold = zero
!
!
!     solve the direct equations (using old code, that we won't touch anymore... )
!     ----------------------------------------------------------------------------
      if ( .not.star ) then
!              
!       loop until max number of iterations
        do it = 1,nitmax
!              
!         1st iteration flag 
          first = ( it.eq.1 )
!          
!         initialize
          vlm(:) = zero
!              
      !$omp parallel do private(basloc,vcos,vsin,vplm) default(shared) &
      !$omp private(isph,pot,vlm,delta) schedule(dynamic,1) 
!              
!         loop over atoms
          do isph = 1,nsph
!              
            call calcv2( first, isph, pot, sigold, basloc, vplm, vcos, vsin )
!            
!           integrate rhs against spherical harmonic Y_l^m
            call intrhs( isph, pot, vlm )
!
!           add potential 
            vlm(:) = vlm(:) + phi(:,isph)
!            
!           solve for sigma^Jac,n
            call solve( isph, vlm, sigma(:,isph) )

!           compute || sigma^n - sigma^n-1 ||_{ S^2 , -1/2 }
            delta(:) = sigma(:,isph) - sigold(:,isph)
            call HSNorm( lmax, delta, norm(isph) )
!            
          end do
!          
      !$omp end parallel do
!          
!         compute rms- and max-norm of "norm" 
          call RMSVec( nsph, norm, drms, dmax )
!
!         decide whether to use DIIS 
          if ( drms.le.tredis .and. ndiis.gt.0 )  dodiis = .true.
!          
!         DIIS solve
          if (dodiis) then
!                  
            xdiis(:,:,nmat) = sigma
            ediis(:,:,nmat) = sigma - sigold
!            
!           iterative solver
            call diis( nbasis*nsph, nmat, xdiis, ediis, bmat, sigma )
!            
          end if
!            
!         printing
          if (iprint.ge.1) then
            ene = pt5*fep*sprod(nbasis*nsph,sigma,psi)
            write(iout,1000) it, ene, drms, dmax
          end if
!
!         desired tolerance has been achieved
          if ( drms.le.tol )  goto 900
!
!         update sigma
          sigold = sigma
!          
        end do
!          
  900   continue
!          
!       compute energy = < psi, sigma >
        ene = sprod( nbasis*nsph, sigma, psi )
!
!
!     solve adjoint equations
!     -----------------------
      else
!
!       loop until maximum number of iterations
        do it = 1,nitmax
!
!         1st iteration flag
          first = ( it.eq.1 )
!            
!         initialize rhs vector
          vlm(:) = zero
!
!         loop over domains
          do isph = 1,nsph
!          
            vlm(:) = psi(:,isph)
!      
!           compute rhs
            call AdjRHS( iout, iprint, first, isph, nbasis, lmax, ngrid, w, csph, &
                         rsph, grid, fi, basis, eta, nl, inl, sigold, vlm, basloc, & 
                         vplm, vcos, vsin, fact )
!                         
!           solve
            call solve( isph, vlm, sigma(:,isph) )
!
!           compute || sigma^n - sigma^n-1 ||_{ S^2 , -1/2 }
            delta(:) = sigma(:,isph) - sigold(:,isph)
            call HSNorm( lmax, delta, norm(isph) )
!            
          end do

!         compute rms- and max-norm of "norm"
          call rmsvec(nsph,norm,drms,dmax)
!          
          if ( drms.le.tredis .and. ndiis.gt.0 )  dodiis=.true.
!
!         DIIS solver
          if (dodiis) then
!
!           initial guess
            xdiis(:,:,nmat) = sigma
!
!           error
            ediis(:,:,nmat) = sigma - sigold
!            
!           iterative solver
            call diis( nbasis*nsph, nmat, xdiis, ediis, bmat, sigma )
!            
          end if
!            
!         printing
          if ( iprint.ge.1 ) then
            write(iout,1000) it, zero, drms, dmax
          end if
!
!         desired tolerance has been reached
          if ( drms.le.tol )  goto 910
!
!         update sigma
          sigold = sigma
!
        end do
!        
  910 continue
!  
      end if
!
!     free the memory
      deallocate( pot, vlm, sigold, delta, norm, vplm, basloc, vcos, vsin, &
                  xdiis, ediis, bmat, scr , stat=istatus )
      if ( istatus.ne.0 ) then            
        write(*,*) 'itsolv2 : deallocation failed !'
        stop
      endif
!
 1000 format(' Energy at iteration ',i4,': ',f14.7,' error (rms,max):',2f14.7)
      return
!
!
endsubroutine itsolv2
