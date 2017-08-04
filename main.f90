!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
!
!
program main
!
      use ddcosmo , only : memuse, memmax, args, nsph, read_x, read_y, read_z, &
                           read_r, read_q, read_control_file, read_molecule_file, &
                           reset_ngrid, ddinit, ncav, nbasis, ccav, iscrf, &
                           igrad, tokcal, memfree
!      
      implicit none
!
!     - point charges [ THIS WILL NEED TO BE CLEANED UP ]
!
      real*8, allocatable :: charge(:)
!
!     - electrostatic potential phi(ncav) and psi vector psi(nbasis,nsph)
!
      real*8, allocatable :: phi(:), psi(:,:)
!
!     - ddcosmo solution sigma(nbasis,nsph) and adjoint solution s(nbasis,nsph)
!
      real*8, allocatable :: sigma(:,:), s(:,:)
!
!     - ddpcm solution phi_eps(nbasis,nsph)   
! 
      real*8, allocatable :: phi_eps(:,:)
!
!     - forces ( cosmo : fx(3,nsph) , pcm : fx (nsph,3) ) [ THIS WILL NEED TO BE CLEANED UP ]
!
      real*8, allocatable :: fx(:,:)
!
!     - miscellanea
!
      integer :: i, istatus, idec, iidec, igrid
      real*8  :: esolv, xx(1)
!
      logical :: interactive_mode = .true.
!
!-------------------------------------------------------------------------------
!
!     initialize
      memuse = 0 ; memmax = 0
!
!     check number of variables in script file
      if ( iargc().ne.2 ) stop
!
!     read control variables
      do i = 1,iargc()
        call getarg( i, args(i) )
      enddo
!
!     read control file
      call read_control_file() 
!
!     adjust number of grid points so that 2*lmax is integrated exactly
      call reset_ngrid( igrid )
!
!     read molecule file
      call read_molecule_file()
!
!     redirect read_q to charge [ THIS WILL NEED TO BE CLEANED UP ]
      allocate( charge(nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'main : [1] allocation failed !'
        stop
      endif
      charge = read_q
!
!     initialize
      call ddinit( nsph, read_x, read_y, read_z, read_r )
!      
!
! ==============================  M O D I F Y   H E R E  ============================== 
!
! Place here your favorite routine to assemble the solute's electrostatic potential
! "phi" and the "psi" vector. Such a routine should replace "mkrhs".
! for classical solutes, assembling the psi vector is straightforward; for qm solutes
! it requires a numerical integration (similar to the one used to compute the xc 
! contributions in dft), as detaild in J. Chem. Phys. 141, 184108
! here, we compute the potential and the psi vector using the supplied routine mkrhs,
! which needs to be replaced by your routine.
!
! -------------------------------------------------------------------------------------
!
!     allocate workspaces
      allocate( phi(ncav), psi(nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'main : [2] allocation failed !'
        stop
      endif
!
!     update memory usage
      memuse = memuse + ncav + nbasis*nsph
      memmax = max(memmax,memuse)
!
!     compute electrostatic potential phi and psi vector
      call mkrhs( nsph, charge, read_x, read_y, read_z, ncav, ccav, phi, nbasis, psi )
!
! ===============================  E N D   M O D I F Y  =============================== 
!
!
!
!======================================================================================
!  I N T E R A C T I V E   M O D E                                                    |
!======================================================================================
!
      if ( interactive_mode ) then
!
!       display option menu in infinite loop 
        do
!         
          write(*,*)'========================================'
          write(*,*)'COSMO .................................1'
          write(*,*)'COSMO & FORCES ........................2'
          write(*,*)'PCM & FORCES ..........................3'
          write(*,*)''
          write(*,*)'DEBUG TESTS ..........................-1'
          write(*,*)'QUIT ..................................0'
          write(*,*)'========================================'
!
          read(*,*) idec
!
!         menu selection
          select case(idec)
!
!         debug tests
!         ===========
          case(-1)
!
            call debug_tests()
!
!         quit
!         ====
          case(0)
!
!           break free of infinite loop
            exit
!
!         cosmo
!         =====  
          case(1)
!                  
!           allocate workspace
            allocate( sigma(nbasis,nsph) , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : COSMO allocation failed !'
              stop
            endif
!
!           solve cosmo equation
            call cosmo( .false., .true., phi, xx, psi, sigma, esolv )
            write (6,'(1x,a,f14.6)') 'ddcosmo electrostatic solvation energy (kcal/mol):', esolv*tokcal
!            
!           deallocate workspace
            deallocate( sigma , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : COSMO deallocation failed !'
              stop
            endif 
!
!         cosmo & forces
!         ==============
          case(2)
!                  
!           allocate workspaces
            allocate( sigma(nbasis,nsph), s(nbasis,nsph), fx(3,nsph) , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : COSMO & FORCES allocation failed !'
              stop
            endif 
!
!           solve cosmo equation
            call cosmo( .false., .true., phi, xx, psi, sigma, esolv )
            write (6,'(1x,a,f14.6)') 'ddcosmo electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
!           solve cosmo adjoint equation [ esolv is not touchedi! ]
            call cosmo( .true., .false., phi, xx, psi, s, esolv )
!
!           compute forces
            call forces( nsph, charge, phi, sigma, s, fx )
!
!           check forces
            write(*,*)'check forces ? 1 - Yes'
            read(*,*) iidec
            if (iidec.eq.1)  call check_forcesCOSMO( esolv, charge, fx )
!             
!           deallocate workspaces
            deallocate( sigma, s, fx , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : COSMO & FORCES deallocation failed !'
              stop
            endif 
!
!         pcm & forces
!         ============
          case(3)
!                  
!           allocate workspaces
            allocate( sigma(nbasis,nsph), phi_eps(nbasis,nsph), fx(nsph,3) , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : PCM & FORCES allocation failed !'
              stop
            endif 
!                  
!           solve pcm equations
            call pcm(   .false.,  .true., .true., phi,      xx, phi_eps )
            call cosmo( .false., .false.,          xx, phi_eps, psi, sigma, esolv )
            write (6,'(1x,a,f14.6)') 'ddpcm electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
!           compute forces
            call compute_forces( phi, charge, psi, sigma, phi_eps, fx )
!
!           check forces
            write(*,*)'check forces ? 1 - Yes'
            read(*,*) iidec
            if (iidec.eq.1)  call check_forcesPCM( charge, fx, esolv )
!            
!           deallocate workspaces
            deallocate( sigma, phi_eps, fx , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : PCM & FORCES deallocation failed !'
              stop
            endif 
!
          endselect
        enddo
!
!
!======================================================================================
!  N O N - I N T E R A C T I V E   M O D E                                            |
!======================================================================================
!
      else
!
!
!       allocate solution vector sigma
        allocate( sigma(nbasis,nsph) , stat=istatus )
        if ( istatus.ne.0 ) then
          write(*,*)'main : [3] allocation failed !'
          stop
        endif
!
!       update memory usage
        memuse = memuse + nbasis*nsph
        memmax = max(memmax,memuse)
!
!
!       COSMO
!       =====
        if ( iscrf.eq.0 ) then
!                
!         solve cosmo equation
          call cosmo( .false., .true., phi, xx, psi, sigma, esolv )
          write (6,'(1x,a,f14.6)') 'ddcosmo electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
!         forces
          if ( igrad.eq.1 ) then
!            
!           allocate workspaces
            allocate( s(nbasis,nsph), fx(3,nsph) , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : [4] allocation failed !'
              stop
            endif 
!            
!           update memory usage
            memuse = memuse + nbasis*nsph + 3*nsph
            memmax = max(memmax,memuse)
!
!           solve cosmo adjoint equation
            call cosmo( .true., .false., phi, xx, psi, s, esolv )
!
!           compute forces
            call forces( nsph, charge, phi, sigma, s, fx )
!!!            call check_forcesCOSMO( esolv, charge, fx )
!             
!           deallocate workspaces
            deallocate( s, fx , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'main : [1] deallocation failed !'
              stop
            endif 
!            
!           update memory usage
            memuse = memuse - nbasis*nsph - 3*nsph
            memmax = max(memmax,memuse)
!
          endif
!
!
!       PCM        
!       ===
        else
!
!         allocate workspaces
          allocate( phi_eps(nbasis,nsph), fx(nsph,3) , stat=istatus )
          if ( istatus.ne.0 ) then
            write(*,*)'main : [5] failed allocation !'
            stop
          endif 
!          
!         update memory usage
          memuse = memuse + nbasis*nsph + 3*nsph
          memmax = max(memmax,memuse)
!                
!         solve pcm equations
          call pcm(   .false.,  .true., .true., phi,      xx, phi_eps )
          call cosmo( .false., .false.,          xx, phi_eps, psi, sigma, esolv )
          write (6,'(1x,a,f14.6)') 'ddpcm electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
!         compute forces
          call compute_forces( phi, charge, psi, sigma, phi_eps, fx )
!!!          call check_forcesPCM( charge, fx, esolv )
!          
!         deallocate workspaces
          deallocate( phi_eps, fx , stat=istatus )
          if ( istatus.ne.0 ) then
            write(*,*)'main : [2] deallocation failed !'
            stop
          endif 
!          
!         update memory usage
          memuse = memuse - nbasis*nsph - 3*nsph
          memmax = max(memmax,memuse)
!
        endif
!
!       deallocate solution vector sigma
        deallocate( sigma , stat=istatus )
        if ( istatus.ne.0 ) then
          write(*,*)'main : [3] deallocation failed !'
          stop
        endif 
!
!     end of non-interactive mode
      endif
!======================================================================================
!
!
!     deallocate workspaces
      deallocate( charge, phi, psi , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'main : [4] deallocation failed !'
        stop
      endif 
!      
!     update memory usage
      memuse = memuse - 5*nsph - ncav - 2*nsph*nbasis
      memmax = max(memmax,memuse)
!
!     free memory
      call memfree
!
      write(6,*) 'maximum quantity of memory allocated:', memmax
!
!
endprogram main
