program main
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
! Written by Filippo Lipparini, October 2015
!            Paolo Gatto,       December 2017
!            Filippo Lipparini, March 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
! Sample interface with the ddcosmo module.                                    !
! This program reads the various ddcosmo parameters and the geometrical        !
! information of the molecule from a text file. The solute is described        !
! by point charges, which are also read from a file.                           !
! Finally, the van der Waals radii are also read in.                           !
!                                                                              !
! This program shows how to call the ddcosmo routines and how to initialize    !
! the various ddcosmo quantities; also, shows where the user should place      !
! his/her routines to compute the molecular electrostatic quantities.          !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
!
integer :: i, ii, isph, ig, n
real*8  :: tobohr, esolv, xx(1)
real*8, parameter :: toang=0.52917721092d0, tokcal=627.509469d0
!
! quantities to be allocated by the user.
! - solute's parameters, such as coordinates, vdw radii and
!   point charges used to model the solute (or multipoles, or
!   qm density...)
!
real*8, allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)
!
! - electrostatic potential phi(ncav) and psi vector psi(nylm,n)
!
real*8, allocatable :: phi(:), psi(:,:)
!
! - ddcosmo solution sigma (nylm,n) and adjoint solution s(nylm,n)
!
real*8, allocatable :: sigma(:,:), s(:,:)
!
! - forces:
!
real*8, allocatable :: fx(:,:), zeta(:), ef(:,:)
!
! - for qm solutes, fock matrix contribution.
!
! here, we read all the ddcosmo parameters from a file named Input.txt
!
open (unit=100,file='Input.txt',form='formatted',access='sequential')
!
! scalar parameters. the variables are defined in the ddcosmo module and are common to
! all the ddcosmo routines (no need to declare them if ddcosmo.mod is loaded.)
!
read(100,*) iprint      ! printing flag
read(100,*) nproc       ! number of openmp threads
read(100,*) lmax        ! max angular momentum of spherical harmonics basis
read(100,*) ngrid       ! number of lebedev points
read(100,*) iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
read(100,*) igrad       ! whether to compute (1) or not (0) forces
read(100,*) eps         ! dielectric constant of the solvent
read(100,*) eta         ! regularization parameter
!
read(100,*) n           ! number of atoms
!
allocate (x(n),y(n),z(n),rvdw(n),charge(n))
!
! we also read from the same file the charges, coordinates and vdw radii.
! in this example, the coordinates and radii are read in angstrom and
! converted in bohr before calling ddinit.
!
do i = 1, n
  read(100,*) charge(i), x(i), y(i), z(i), rvdw(i)
end do
tobohr = 1.0d0/toang
x    = x*tobohr
y    = y*tobohr
z    = z*tobohr
rvdw = rvdw*tobohr
!
close (100)
!
! call the initialization routine. this routine allocates memory, computes some
! quantities for internal use and creates and fills an array ccav(3,ncav) with
! the coordinates of the grid points at which the user needs to compute the potential.
! ncav is the number of external grid points and nylm the number of spherical
! harmonics functions used for the expansion of the various ddcosmo quantities;
! both are computed by ddinit and defined as common variables in ddcosmo.mod.
!
call ddinit(n,x,y,z,rvdw)
!
allocate (phi(ncav),psi(nylm,n))
!
! --------------------------   modify here  --------------------------  
!
! place here your favorite routine to assemble the solute's electrostatic potential
! and the "psi" vector. Such a routine should replace "mkrhs".
! for classical solutes, assembling the psi vector is straightforward; for qm solutes
! it requires a numerical integration (similar to the one used to compute the xc 
! contributions in dft), as detaild in J. Chem. Phys. 141, 184108
! here, we compute the potential and the psi vector using the supplied routine mkrhs,
! which needs to be replaced by your routine.
!
call mkrhs(n,charge,x,y,z,ncav,ccav,phi,nylm,psi)
!
! --------------------------   end modify   --------------------------  
!
! now, call the ddcosmo solver
!
allocate (sigma(nylm,n))
!
call cosmo(.false., .true., phi, xx, psi, sigma, esolv)
!
if (iprint.ge.3) call prtsph('solution to the ddCOSMO equation',nsph,0,sigma)
!
write (6,'(1x,a,f14.6)') 'ddcosmo electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
! this is all for the energy. if the forces are also required, call the solver for
! the adjoint problem. 
! the solution to the adjoint system is required also to compute the Fock matrix 
! contributions.
!
if (igrad.eq.1) then
  write(6,*)
  allocate (s(nylm,n))
  allocate (fx(3,n))
  call cosmo(.true., .false., xx, xx, psi, s, esolv)
!
  if (iprint.ge.3) call prtsph('solution to the ddCOSMO adjoint equation',nsph,0,s)
!
! now call the routine that computes the ddcosmo specific contributions to the forces.
!
  call forces_dd(n,phi,sigma,s,fx)
!
! form the "zeta" intermediate
!
  allocate (zeta(ncav))
  call ddmkzeta(s,zeta)
!
  if (iprint.ge.4) call ptcart('zeta',nsph,0,zeta)
!
! --------------------------   modify here  --------------------------  
!
! finally, add the contributions that depend on the derivatives of the potential
! on the derivatives of the potential, i.e., the electric field
! produced by the solute at the cavity points times the ddcosmo intermediate zeta 
! and the solute's potential derivatives at the cavity points times zeta. 
! the two terms are described in JCP, 141, 184108, eqs. 47, 48.
!
! for a solute represented by point charges, the two terms can be rearranged as zeta 
! contracted with the electric field of the solute plus the electric field produced
! by zeta, which have the physical dimension of a charge, at the nuclei times the 
! charges. This rearrangement allows one to use fast summations methods, such as the
! fast multipole method, for this task.
!
! a routine for point charges that follows the aformentioned strategy is provided as
! an example (efld). note that the coordinates should be packed into an array of
! dimension (3,n). we use csph, as it contains exactly this.
! the user will need to replace efld with his/her favorite routine.
!
  allocate(ef(3,max(n,ncav)))
!
! 1. solute's electric field at the cav points times zeta:
!
!   compute the electric field
!
  call efld(n,charge,csph,ncav,ccav,ef)
!
!   contract it with the zeta intermediate
!
  ii = 0
  do isph = 1, nsph
    do ig = 1, ngrid
      if (ui(ig,isph).gt.zero) then 
        ii = ii + 1
        fx(:,isph) = fx(:,isph) - zeta(ii)*ef(:,ii)
      end if
    end do
  end do
!
! 2. "zeta's" electric field at the nuclei times the charges. 
!
!   compute the "electric field"
!
  call efld(ncav,zeta,ccav,n,csph,ef)
!
!   contract it with the solute's charges.
!
  do isph = 1, nsph
    fx(:,isph) = fx(:,isph) - ef(:,isph)*charge(isph)
  end do
!
! for point charges, there is no contribution from the derivatives of the psi vector.
! for quantum mechanical solutes, such a contribution needs to be handled via a numerical
! integration.
!
! --------------------------   end modify   --------------------------  
! 
  deallocate (zeta, ef)
!
  if (iprint.ge.1) then
    write(iout,2000)
2000 format(1x,'ddCOSMO forces (atomic units):',/, &
              1x,' atom',15x,'x',15x,'y',15x,'z')
    do isph = 1, nsph
      write(6,'(1x,i5,3f16.8)') isph, fx(:,isph)
    end do
  end if
end if
!
! clean up:
!
deallocate (x,y,z,rvdw,charge,phi,psi,sigma)
!
if (igrad.eq.1) deallocate (s,fx)
call memfree
!
end program main
