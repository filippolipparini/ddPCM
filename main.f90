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
! grids by D. Laikov and C. van Wuellen, as publicly available on CCL. If the routines
! are used, the following reference should also be included:
!
! [6] V.I. Lebedev, and D.nsph. Laikov
!     "A quadrature formula for the sphere of the 131st
!      algebraic order of accuracy"
!     Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
! Written by Filippo Lipparini, October 2015.
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
integer :: i, istatus
real*8  :: tobohr, esolv
real*8, parameter :: toang=0.52917721092d0, tokcal=627.509469d0
!
! quantities to be allocated by the user.
! - solute's parameters, such as coordinates, vdw radii and
!   point charges used to model the solute (or multipoles, or
!   qm density...)
!
real*8, allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)
!
! - electrostatic potential phi(ncav) and psi vector psi(nbasis,nsph)
!
real*8, allocatable :: phi(:), psi(:,:), g(:,:)
!
! - ddcosmo solution sigma (nbasis,nsph) and adjoint solution s(nbasis,nsph)
!
real*8, allocatable :: sigma(:,:), s(:,:)
!
! - forces:
!
real*8, allocatable :: fx(:,:)
!
! - for qm solutes, fock matrix contribution.
!
character(len=64), dimension(2) :: args
!
!-----------------------------------------------------------
! here, we read all the ddcosmo parameters from a file named Input.txt
!
memuse = 0
memmax = 0
!
      if ( iargc().ne.2 ) stop
      do i = 1,iargc()
        call getarg( i, args(i) )
      enddo
!
!     read control parameters
!     -----------------------
      open( unit=10, file=args(1) )
      read(10,*) iprint      ! printing flag
      read(10,*) nproc       ! number of openmp threads
      read(10,*) lmax        ! max angular momentum of spherical harmonics basis
      read(10,*) ngrid       ! number of lebedev points
      read(10,*) iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
      read(10,*) igrad       ! whether to compute (1) or not (0) forces
      read(10,*) iscrf       ! whether to use cosmo (0) or pcm (1)
      read(10,*) eps         ! dielectric constant of the solvent
      read(10,*) iunit       ! whether to convert to bohr (0) or not (1)
      read(10,*) eta, se     ! regularization parameter
      close(10)
!
      call reset_ngrid
!
!     read atoms file
!     ---------------
      open( unit=10, file=args(2) )
      read(10,*) nsph           ! number of atoms
!
      allocate( x(nsph), y(nsph), z(nsph), rvdw(nsph), charge(nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'main : [1] failed allocation !'
        stop
      endif
!
      memuse = memuse + 5*nsph
      memmax = max(memmax,memuse)
!
      do i = 1, nsph
        read(10,*) charge(i), x(i), y(i), z(i), rvdw(i)
      end do
      tobohr = 1.0d0
      if (iunit.eq.0) tobohr = 1.0d0/toang
      x    = x*tobohr
      y    = y*tobohr
      z    = z*tobohr
      rvdw = rvdw*tobohr
!
      close (10)
!
! call the initialization routine. this routine allocates memory, computes some
! quantities for internal use and creates and fills an array ccav(3,ncav) with
! the coordinates of the grid points at which the user needs to compute the potential.
! ncav is the number of external grid points and nbasis the number of spherical
! harmonics functions used for the expansion of the various ddcosmo quantities;
! both are computed by ddinit and defined as common variables in ddcosmo.mod.
!
      call ddinit(nsph,x,y,z,rvdw)
      allocate( phi(ncav), psi(nbasis,nsph), g(ngrid,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'main : [2] failed allocation !'
        stop
      endif
!
      memuse = memuse + ncav + nbasis*nsph
      memmax = max(memmax,memuse)
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
      call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
!
! --------------------------   end modify   --------------------------  
!
! now, call the ddcosmo solver
!
      allocate( sigma(nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'main : [3] failed allocation !'
        stop
      endif

!
      memuse = memuse + nbasis*nsph
      memmax = max(memmax,memuse)
!
!     COSMO
!     =====
      if (iscrf.eq.0) then
!1              
!       compute \sigma and solvation energy
        call itsolv( .false., phi, psi, sigma, esolv )
        write (6,'(1x,a,f14.6)') 'ddcosmo electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
!     PCM        
!     ===
      else
!              
        call wghpot(phi,g)
        call iefpcm(g,psi,sigma)
!        
      end if
!
! this is all for the energy. if the forces are also required, call the solver for
! the adjoint problem. 
! the solution to the adjoint system is required also to compute the Fock matrix 
! contributions.
!
if (igrad.eq.1) then
  write(6,*)
  allocate( s(nbasis,nsph), fx(3,nsph) , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*)'main : [4] failed allocation !'
    stop
  endif 
  memuse = memuse + nbasis*nsph + 3*nsph
  memmax = max(memmax,memuse)
  call itsolv(.true.,phi,psi,s,esolv)
!
! now call the routine that computes the forces. such a routine requires the potential 
! derivatives at the cavity points and the electric field at the cavity points: it has
! therefore to be personalized by the user. it is included in this sample program as
! forces.f90.
!
   call forces(nsph,charge,phi,sigma,s,fx)
end if
!
! clean up:
!
deallocate( x, y, z, rvdw, charge, phi, psi, sigma , stat=istatus)
if ( istatus.ne.0 ) then
  write(*,*)'main : [1] failed deallocation !'
  stop
endif 

memuse = memuse - 5*nsph - ncav - 2*nsph*nbasis
memmax = max(memmax,memuse)
!
if (igrad.eq.1) then
  deallocate( s, fx , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*)'main : [2] failed deallocation !'
    stop
  endif 

endif
call memfree
!
write(6,*) 'maximum quantity of memory allocated:', memmax
!
end program main
