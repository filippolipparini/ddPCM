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
!
!

subroutine lx(n, x, y)
use ddcosmo
implicit none 
!
! given a vector x, compute y = Lx, where L is the ddCOSMO matrix.
! if dodiag is set to .true., L includes the diagonal blocks, otherwise
! L only includes the off-diagonal ones.
!
integer, intent(in) :: n
real*8,  dimension(nbasis,nsph), intent(in)    :: x
real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
integer             :: isph, istatus
real*8, allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
!
! allocate some memory:
!
allocate (pot(ngrid), vplm(nbasis), basloc(nbasis), vcos(lmax+1), vsin(lmax+1), stat=istatus)
if (istatus .ne. 0) then
  write(*,*) ' lx: allocation failed.'
end if
!
y = zero
!
!$omp parallel do default(shared) private(isph,pot,basloc,vplm,vcos,vsin) &
!$omp schedule(dynamic)
do isph = 1, nsph
  call calcv2(.false., isph, pot, x, basloc, vplm, vcos, vsin)
  call intrhs(isph, pot, y(:,isph))
  y(:,isph) = - y(:,isph)
  if (do_diag) y(:,isph) = y(:,isph) + x(:,isph)/facl
end do
!
deallocate (pot, basloc, vplm, vcos, vsin)
!
return
end
!
subroutine lstarx(n, x, y)
use ddcosmo
implicit none 
!
! given a vector x, compute y = L*x, where L* is the adjoint ddCOSMO matrix.
! if dodiag is set to .true., L includes the diagonal blocks, otherwise
! L only includes the off-diagonal ones.
!
integer, intent(in) :: n
real*8,  dimension(nbasis,nsph), intent(in)    :: x
real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
integer             :: isph, ig, istatus
real*8, allocatable :: xi(:,:), vplm(:), basloc(:), vcos(:), vsin(:)
!
! allocate some memory:
!
allocate (xi(ngrid,nsph), vplm(nbasis), basloc(nbasis), vcos(lmax+1), vsin(lmax+1), stat=istatus)
if (istatus .ne. 0) then
  write(*,*) ' lx: allocation failed.'
end if
!
y = zero
!
!$omp parallel do default(shared) private(isph,ig)
do isph = 1, nsph
  do ig = 1, ngrid
    xi(ig,isph) = dot_product(x(:,isph), basis(:,ig))
  end do
end do
!
!$omp parallel do default(shared) private(isph,basloc,vplm,vcos,vsin) &
!$omp schedule(dynamic)
do isph = 1, nsph
  call adjrhs1(isph, xi, y(:,isph), basloc, vplm, vcos, vsin)
  y(:,isph) = - y(:,isph)
  if (do_diag) y(:,isph) = y(:,isph) + x(:,isph)/facl
end do
!
deallocate (xi, basloc, vplm, vcos, vsin)
!
return
end
!
subroutine ldm1x(n, x, y)
use ddcosmo
implicit none
!
! given a vector x, apply the inverse diagonal of the L matrix:
!
integer,                         intent(in)    :: n
real*8,  dimension(nbasis,nsph), intent(in)    :: x
real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
integer                                        :: isph
!
do isph = 1, nsph
  y(:,isph) = facl*x(:,isph)
end do
!
return
end
!
real*8 function hnorm(n,x)
use ddcosmo
implicit none
!
! compute the h^-1/2 norm of the increment on each sphere, then take the
! rms value.
!
integer,                         intent(in) :: n
real*8,  dimension(nbasis,nsph), intent(in) :: x
!
integer                                     :: isph, istatus
real*8                                      :: vrms, vmax
real*8, allocatable                         :: u(:)
!
allocate (u(nsph),stat=istatus)
if (istatus .ne. 0) then
  write(*,*) ' hnorm: allocation failed.'
  stop
end if
!
do isph = 1, nsph
  call hsnorm(x(:,isph),u(isph))
end do
!
call rmsvec(nsph,u,vrms,vmax)
!
hnorm = vrms
return
!
end
!
subroutine plx(n, x, y)
use ddcosmo
implicit none 
!
! given a vector x, compute y = Lx, where L is the ddCOSMO matrix, then
! apply the inverse diagonal as a preconditioner.
!
integer, intent(in) :: n
real*8,  dimension(nbasis,nsph), intent(in)    :: x
real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
do_diag = .true.
call lx(n,x,y)
call ldm1x(n,y,y)
!
return
end
!
subroutine plstarx(n, x, y)
use ddcosmo
implicit none 
!
! given a vector x, compute y = L*x, where L* is the adjoint ddCOSMO matrix,
! then apply the inverse diagonal as a preconditioner.
!
integer, intent(in) :: n
real*8,  dimension(nbasis,nsph), intent(in)    :: x
real*8,  dimension(nbasis,nsph), intent(inout) :: y
!
do_diag = .true.
call lstarx(n,x,y)
call ldm1x(n,y,y)
!
return
end
!
