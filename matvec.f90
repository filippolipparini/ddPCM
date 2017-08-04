! 
!      888      888  .d8888b.   .d88888b.   .d8888b.  888b     d888  .d88888b.  
!      888      888 d88P  Y88b d88P" "Y88b d88P  Y88b 8888b   d8888 d88P" "Y88b 
!      888      888 888    888 888     888 Y88b.      88888b.d88888 888     888 
!  .d88888  .d88888 888        888     888  "Y888b.   888Y88888P888 888     888 
! d88" 888 d88" 888 888        888     888     "Y88b. 888 Y888P 888 888     888 
! 888  888 888  888 888    888 888     888       "888 888  Y8P  888 888     888 
! Y88b 888 Y88b 888 Y88b  d88P Y88b. .d88P Y88b  d88P 888   "   888 Y88b. .d88P 
!  "Y88888  "Y88888  "Y8888P"   "Y88888P"   "Y8888P"  888       888  "Y88888P"  
!                                                                              
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  COPYRIGHT (C) 2015 by Filippo Lipparini, Benjamin Stamm, Eric Cancès,       !
!  Yvon Maday, Jean-Philip Piquemal, Louis Lagardère and Benedetta Mennucci.   !
!                             ALL RIGHT RESERVED.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! A modular implementation of COSMO using a domain decomposition linear scaling
! strategy.
!
! This code is governed by the LGPL license and abiding by the rules of 
! distribution of free software.
! This program is distributed in the hope that it will be useful, but  
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Lesser General Public License for more details.
!
! Users of this code are asked to include the following references in their
! publications:
!
! [1] E. Cancès, Y. Maday, B. Stamm
!     "Domain decomposition for implicit solvation models"
!     J. Chem. Phys. 139, 054111 (2013)
!
! [2] F. Lipparini, B. Stamm, E. Cancès, Y. Maday, B. Mennucci
!     "Fast Domain Decomposition Algorithm for Continuum Solvation Models: 
!      Energy and First Derivatives"
!     J. Chem. Theory Comput. 9, 3637–3648 (2013)
!
! Also, include one of the three following reference depending on whether you
! use this code in conjunction with a QM [3], Semiempirical [4] or Classical [5]
! description of the solute:
!
! [3] F. Lipparini, G. Scalmani, L. Lagardère, B. Stamm, E. Cancès, Y. Maday,
!     J.-P. Piquemal, M. J. Frisch, B. Mennucci
!     "Quantum, classical, and hybrid QM/MM calculations in solution: General 
!      implementation of the ddCOSMO linear scaling strategy"
!     J. Chem. Phys. 141, 184108 (2014)
!     (for quantum mechanical models)
!
! [4] F. Lipparini, L. Lagardère, G. Scalmani, B. Stamm, E. Cancès, Y. Maday,
!     J.-P. Piquemal, M. J. Frisch, B. Mennucci
!     "Quantum Calculations in Solution for Large to Very Large Molecules: 
!      A New Linear Scaling QM/Continuum Approach"
!     J. Phys. Chem. Lett. 5, 953-958 (2014)
!     (for semiempirical models)
!
! [5] F. Lipparini, L. Lagardère, C. Raynaud, B. Stamm, E. Cancès, B. Mennucci
!     M. Schnieders, P. Ren, Y. Maday, J.-P. Piquemal
!     "Polarizable Molecular Dynamics in a Polarizable Continuum Solvent"
!     J. Chem. Theory Comput. 11, 623-634 (2015)
!     (for classical models, including polarizable force fields
!     
! The users of this code should also include the appropriate reference to the
! COSMO model. This distribution includes the routines to generate lebedev
! grids by D. Laikov and C. van Wuellen, as publicly available on CCL. If the
! routines are used, the following reference should also be included:
!
! [6] V.I. Lebedev, and D.N. Laikov
!     "A quadrature formula for the sphere of the 131st
!      algebraic order of accuracy"
!     Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
! Written by Filippo Lipparini, October 2015.
! Modified by FL and Paolo Gatto, June 2017
!
!-------------------------------------------------------------------------------
! 
! This file collects all the routines that are used to perform a
! matrix-vector multiplication for COSMO and PCM. This includes
!
!   lx      : COSMO matrix
!   lstarx  : COSMO adjoint matrix
!   plx     : COSMO matrix and COSMO inverse diagonal as a preconditioner
!   plstarx : COSMO adjoint matrix and COSMO inverse diagonal as a preconditioner
!
!   rx      : PCM matrix
!   rstarx  : PCM adjoint matrix
!   prx     : PCM matrix and PCM inverse diagonal as a preconditioner
!   prstarx : PCM adjoint matrix and PCM inverse diagonal as a preconditioner
!
!   + service routines
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = Lx, where L is the ddCOSMO matrix.
! if dodiag is set to .true., L includes the diagonal blocks, otherwise
! L only includes the off-diagonal ones.
!-------------------------------------------------------------------------------
!
subroutine lx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, zero, calcv2, intrhs, &
                           do_diag, facl
!      
      implicit none 
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
      !
      integer             :: isph, istatus
      real*8, allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
!      
!-------------------------------------------------------------------------------
!
!     allocate workspaces
      allocate( pot(ngrid), vplm(nbasis), basloc(nbasis), vcos(lmax+1), &
                vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'lx: allocation failed !'
        stop
      endif
!
!     initialize
      y = zero
      !
      !$omp parallel do default(shared) private(isph,pot,basloc,vplm,vcos,vsin) &
      !$omp schedule(dynamic)
      !
!
!     loop over spheres
      do isph = 1,nsph
!      
        call calcv2( .false., isph, pot, x, basloc, vplm, vcos, vsin )
        call intrhs( isph, pot, y(:,isph) )
!
!       P. Gatto : why the sign flip ???
        y(:,isph) = - y(:,isph)
!
!       add action of diagonal block
        if ( do_diag )  y(:,isph) = y(:,isph) + x(:,isph)/facl
!        
      enddo
!
!     deallocate workspaces
      deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'lx: allocation failed !'
        stop
      endif
!
!
endsubroutine lx
!-------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = L*x, where L* is the adjoint ddCOSMO matrix.
! if dodiag is set to .true., L includes the diagonal blocks, otherwise
! L only includes the off-diagonal ones.
!-------------------------------------------------------------------------------
subroutine lstarx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, zero, basis, do_diag, &
                           adjrhs1, facl
!      
      implicit none 
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
      integer             :: isph, ig, istatus
      real*8, allocatable :: xi(:,:), vplm(:), basloc(:), vcos(:), vsin(:)
!
!-------------------------------------------------------------------------------
!
!     allocate workspaces
      allocate( xi(ngrid,nsph), vplm(nbasis), basloc(nbasis), vcos(lmax+1), &
                vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'lstarx: allocation failed!'
        stop
      endif
!
!     initilize
      y = zero
      !
      !$omp parallel do default(shared) private(isph,ig)
      !
!
!     expand x over spherical harmonics
!     ---------------------------------
!
!     loop over spheres      
      do isph = 1,nsph
!
!       loop over gridpoints
        do ig = 1,ngrid
!        
          xi(ig,isph) = dot_product( x(:,isph), basis(:,ig) )
!          
        enddo
      enddo
      !
      !$omp parallel do default(shared) private(isph,basloc,vplm,vcos,vsin) &
      !$omp schedule(dynamic)
!
!     compute action
!     --------------
!
!     loop over spheres
      do isph = 1,nsph
!      
        call adjrhs1( isph, xi, y(:,isph), basloc, vplm, vcos, vsin )
!        
!       P. Gatto : why the sign flip ???
        y(:,isph) = - y(:,isph)
!
!       add action of diagonal block
        if ( do_diag )  y(:,isph) = y(:,isph) + x(:,isph)/facl
!        
      enddo
!
!     deallocate workspaces
      deallocate( xi, basloc, vplm, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'lstarx: allocation failed !'
        stop
      endif
!
!
endsubroutine lstarx
!-------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, apply the inverse diagonal (block) of the L matrix:
!-------------------------------------------------------------------------------
!
subroutine ldm1x( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, facl
!      
      implicit none
!
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
      integer                                        :: isph
!
!-------------------------------------------------------------------------------
!
!     loop over spheres
      do isph = 1,nsph
!
!       apply inverse
        y(:,isph) = facl*x(:,isph)
!        
      enddo
!
!
endsubroutine ldm1x
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
! compute the h^-1/2 norm of the increment on each sphere, then take the
! rms value.
!-------------------------------------------------------------------------------
!
real*8 function hnorm( n, x )
!
      use ddcosmo , only : nbasis, nsph, hsnorm
!
      implicit none
      integer,                         intent(in) :: n
      real*8,  dimension(nbasis,nsph), intent(in) :: x
!
      integer                                     :: isph, istatus
      real*8                                      :: vrms, vmax
      real*8, allocatable                         :: u(:)
!
!-------------------------------------------------------------------------------
!
!     allocate workspace
      allocate( u(nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'hnorm: allocation failed !'
        stop
      endif
!
!     loop over spheres
      do isph = 1,nsph
!
!       compute norm contribution
        call hsnorm( x(:,isph), u(isph) )
      enddo
!
!     compute rms of norms
      call rmsvec( nsph, u, vrms, vmax )
!
!     return value
      hnorm = vrms
!
!     deallocate workspace
      deallocate( u , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'hnorm: deallocation failed !'
        stop
      endif
!
!
endfunction hnorm
!-------------------------------------------------------------------------------
!
!
!
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = Lx, where L is the ddCOSMO matrix, then
! apply the inverse diagonal as a preconditioner.
!-------------------------------------------------------------------------------
!
subroutine plx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, do_diag
!
      implicit none 
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
!-------------------------------------------------------------------------------
!
!     activate action of diagonal blocks
      do_diag = .true.
!      
!     action of COSMO
      call lx( n, x, y )
!
!     action of inverse diagonal
      call ldm1x( n, y, y )
!
!
endsubroutine plx
!-------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = L*x, where L* is the adjoint ddCOSMO matrix,
! then apply the inverse diagonal as a preconditioner.
!-------------------------------------------------------------------------------
!
subroutine plstarx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, do_diag
!
      implicit none 
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
!-------------------------------------------------------------------------------
!
!     activate action of diagonal blocks
      do_diag = .true.
!
!     action of adjoint COSMO
      call lstarx( n, x, y )
!
!     action of inverse diagonal
      call ldm1x( n, y, y )
!
!  
endsubroutine plstarx
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, solve Px = y, where P is a given preconditioner.
! note that we assume that P^{-1} is available in precm1 here...
!-------------------------------------------------------------------------------
!
subroutine precx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, zero, precm1, one
!      
      implicit none
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
      integer :: isph
!
!-------------------------------------------------------------------------------
!
!     loop over spheres
      do isph = 1,nsph
!
!       mutiply by preconditioner
        call DGEMV( 'N', nbasis, nbasis, one, precm1(:,:,isph), nbasis, &
                     x(:,isph), 1, zero, y(:,isph), 1 )
!                     
      enddo
!
!
endsubroutine precx
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = R_\eps x
!-------------------------------------------------------------------------------
!
subroutine rx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, eps
!
      implicit none
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
      integer                           :: isph, istatus
      real*8, allocatable, dimension(:) :: ulm, u, basloc, vplm, vcos, vsin
!
!-------------------------------------------------------------------------------
!
!     allocate workspaces
      allocate( ulm(nbasis), u(ngrid), basloc(nbasis), vplm(nbasis), &
                vcos(lmax+1), vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'rx: allocation failed!'
        stop
      endif
!
!     loop over spheres
      do isph = 1, nsph
!
!       action of PCM
        call mkrvec( isph, eps, x, y(:,isph), ulm, u, basloc, vplm, vcos, vsin )
      enddo
!
!     deallocate workspaces
      deallocate( ulm, u, basloc, vplm, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'rx: deallocation failed!'
        stop
      endif
!
!
endsubroutine rx
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute z = R_\eps x. 
! then, solve y = Py.
!-------------------------------------------------------------------------------
!
subroutine prx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, do_diag
!
      implicit none
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
      integer                              :: istatus
      real*8,  allocatable, dimension(:,:) :: z(:,:)
!      
!-------------------------------------------------------------------------------
!
!     allocate workspace
      allocate( z(nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'prx: allocation failed!'
        stop
      endif
!
!     active action of diagonal block
      do_diag = .true.
!
!     action of PCM
      call rx( n, x, z )
!
!     apply preconditioner
      call precx( n, z, y )
!
!     deallocate workspace
      deallocate( z , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'prx: deallocation failed!'
        stop
      endif
!
!
endsubroutine prx
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = R_\eps^* x
!-------------------------------------------------------------------------------
!
subroutine rstarx( n, x, y )
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, eps
!      
      implicit none
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
      integer                           :: isph, istatus
      real*8, allocatable, dimension(:) :: ulm, u, basloc, vplm, vcos, vsin
!
!-------------------------------------------------------------------------------
!
!     allocate workspaces
      allocate( ulm(nbasis), u(ngrid), basloc(nbasis), vplm(nbasis), &
                vcos(lmax+1), vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'rstarx: allocation failed!'
        stop
      endif
!
!     loop over spheres
      do isph = 1,nsph
!
!       action of adjoint PCM
        call adjvec( isph, eps, x, y(:,isph), ulm, u, basloc, vplm, vcos, vsin )
!        
      enddo
!
!     deallocate workspaces
      deallocate( ulm, u, basloc, vplm, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'rstarx: deallocation failed!'
        stop
      endif
!
!
endsubroutine rstarx
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute z = R_\eps x. 
! then, solve y = Py.
!-------------------------------------------------------------------------------
!
subroutine prstarx(n,x,y)
!
      use ddcosmo , only : nbasis, nsph, do_diag
!      
      implicit none
      integer,                         intent(in)    :: n
      real*8,  dimension(nbasis,nsph), intent(in)    :: x
      real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
      integer                              :: istatus
      real*8,  allocatable, dimension(:,:) :: z(:,:)
!
!-------------------------------------------------------------------------------
!
!     allocate workspace
      allocate( z(nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'prstarx: allocation failed!'
      endif
!
!     activate action of diagonal block
      do_diag = .true.
!
!     action of adjoint PCM
      call rstarx( n, x, z )
!
!     apply preconditioner
      call precx( n, z, y )
!
!     deallocate workspace
      deallocate( z , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'prstarx: deallocation failed!'
      endif
!
!
endsubroutine prstarx
