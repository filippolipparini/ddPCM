subroutine forces_dd(n,phi,sigma,s,fx)
use ddcosmo
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
!  COPYRIGHT (C) 2015 by Filippo Lipparini, Benjamin Stamm, Paolo Gatto        !
!  Eric Cancès, Yvon Maday, Jean-Philip Piquemal, Louis Lagardère and          !
!  Benedetta Mennucci.                                                         !
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
! grids by D. Laikov and C. van Wuellen, as publicly available on CCL. If the routines
! are used, the following reference should also be included:
!
! [6] V.I. Lebedev, and D.N. Laikov
!     "A quadrature formula for the sphere of the 131st
!      algebraic order of accuracy"
!     Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
! Written by Filippo Lipparini, October 2015.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                               !
! Sample driver for the calculation of the ddCOSMO forces.                     !
                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
!
integer,                         intent(in)    :: n
real*8,  dimension(ncav),        intent(in)    :: phi
real*8,  dimension(nylm,nsph),   intent(in)    :: sigma, s
real*8,  dimension(3,n),         intent(inout) :: fx
!
integer :: isph, ig, ii, c1, c2, cr
real*8  :: fep
!
real*8, allocatable :: xi(:,:), phiexp(:,:), zeta(:), ef(:,:)
real*8, allocatable :: basloc(:), dbsloc(:,:), vplm(:), vcos(:), vsin(:)
!
allocate (xi(ngrid,nsph),phiexp(ngrid,nsph))
allocate (basloc(nylm),dbsloc(3,nylm),vplm(nylm),vcos(lmax+1),vsin(lmax+1))
!
! initialize the timer:
!
call system_clock(count_rate=cr)
call system_clock(count=c1)
!
! compute xi:
!
!$omp parallel do default(shared) private(isph,ig)
do isph = 1, nsph
  do ig = 1, ngrid
    xi(ig,isph) = dot_product(s(:,isph),basis(:,ig))
  end do
end do
!$omp end parallel do
!
if (iprint.ge.4) call ptcart('xi',nsph,0,xi)
!
! expand the potential on a sphere-by-sphere basis (needed for parallelism):
!
ii = 0
phiexp = zero
do isph = 1, nsph
  do ig = 1, ngrid
    if (ui(ig,isph).gt.zero) then
      ii = ii + 1
      phiexp(ig,isph) = phi(ii)
    end if
  end do
end do
!
fx = zero
do isph = 1, nsph
  call fdoka(isph,sigma,xi(:,isph),basloc,dbsloc,vplm,vcos,vsin,fx(:,isph)) 
  call fdokb(isph,sigma,xi,basloc,dbsloc,vplm,vcos,vsin,fx(:,isph)) 
  call fdoga(isph,xi,phiexp,fx(:,isph)) 
end do
!
2000 format(1x,'ddCOSMO-only contributions to the forces (atomic units):',/, &
              1x,' atom',15x,'x',15x,'y',15x,'z')
!
if (iprint.ge.4) then
  write(iout,2000)
  do isph = 1, nsph
    write(6,'(1x,i5,3f16.8)') isph, fx(:,isph)
  end do
end if
!
deallocate (basloc,dbsloc,vplm,vcos,vsin)
!
call system_clock(count=c2)
if (iprint.gt.0) then
  write(iout,1010) dble(c2-c1)/dble(cr)
1010 format(' the computation of the ddCOSMO part of the forces took ',f8.3,' seconds.')
end if
!
deallocate (xi,phiexp)
!
! scale the forces time the cosmo factor:
!
fep = pt5*(eps-one)/eps
fx  = fep*fx
!
return
end
