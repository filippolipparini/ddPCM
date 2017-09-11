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
! [6] V.I. Lebedev, and D.N. Laikov
!     "A quadrature formula for the sphere of the 131st
!      algebraic order of accuracy"
!     Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
! Written by Filippo Lipparini, October 2015.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
! Sample driver for the calculation of the ddCOSMO forces.                     !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------
! Recall :
!
!  E_s = 1/2 f(eps) < Psi , sigma >  ;  F_i = - dE_s / dr_i  ;  L^T s = Psi
!
! Then, if we indicate a generic derivative with ', we obtain :
!
!   g' = (L sigma)' = L' sigma + L sigma'
!
! Thus :
!
!   E_s' = ... <   Psi , sigma' >
!
!        = ... < L^T s , sigma' >
!
!        = ... <     s , L sigma' >
!
!        = ... < s , g' > - < s , L' sigma >
!
! Hence :
!
!   F = ... ( - < s , g' > + < s , L' sigma > )
!
! Recall :
!
!   < s , L' sigma > = K_a + K_b  [ see PCM notes ]
!
! and
!
!   < s , g' > = sum sum sum  Y_l^m(y_n) [s_j]_l^m  w_n   ( U_n^j' Phi_n^j + U_n^j  Phi_n^j' )
!                 j   n  l,m
!
!                        |------- xi(j,n) -------|
!
!
!              = sum sum  xi(j,n) w_n U_n^j  Phi_n^j'  +  xi(j,n) w_n U_n^j' Phi_n^j
!                 j   n
!
!                         |-- zeta(j,n) --|
!
! Finally :
!
!   Phi_n^i' = \grad_i Phi_n^i = -E[z_i]_n^i
!
! where E... is the electric field produced by all charges z_j such that
! j \ne i , and  
!
!   Phi_n^j' = \grad_i Phi_n^j = -E(z_i)_n^j 
!
! where E... is the electric field produced by the charge z_i at the cavity
! points for j \ne i .
!-------------------------------------------------------------------------------
!
subroutine forces( n, charge, phi, sigma, s, fx )
!
      use ddcosmo
!
      implicit none
!
      integer,                         intent(in)    :: n
      real*8,  dimension(n),           intent(in)    :: charge
      real*8,  dimension(ncav),        intent(in)    :: phi
      real*8,  dimension(nbasis,nsph), intent(in)    :: sigma
      real*8,  dimension(nbasis,nsph), intent(in)    :: s
      real*8,  dimension(3,n),         intent(inout) :: fx
!
      integer :: isph, ig, ii, c1, c2, cr, istatus
      real*8  :: fep
!
      real*8, allocatable :: xi(:,:), phiexp(:,:), zeta(:), ef(:,:)
      real*8, allocatable :: basloc(:), dbsloc(:,:), vplm(:), vcos(:), vsin(:)
!      
!-------------------------------------------------------------------------------
!
!     allocate workspaces...
      allocate( xi(ngrid,nsph), phiexp(ngrid,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'forces : [1] failed allocation !'
      endif
!
!     ... and more workspaces      
      allocate( basloc(nbasis), dbsloc(3,nbasis), vplm(nbasis), vcos(lmax+1), & 
                vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'forces : [2] failed allocation !'
      endif
!
!     initialize the timer
      call system_clock( count_rate = cr )
      call system_clock( count = c1 )
!
      !$omp parallel do default(shared) private(isph,ig)
!
!     loop over atoms
      do isph = 1, nsph
!
!       loop over gridpoints
        do ig = 1, ngrid
!        
!         compute xi(j,n)
!         ===============
          xi(ig,isph) = dot_product( s(:,isph), basis(:,ig) )
!          
        enddo
      enddo
!      
      !$omp end parallel do
!
!
!     initialize
      ii = 0 ; phiexp = zero
!
!     loop over atoms
      do isph = 1, nsph
!
!       loop over integration points
        do ig = 1, ngrid
!
!         positive contribution from integration point
          if ( ui(ig,isph).gt.zero ) then
!                  
!           advance index
            ii = ii + 1
!
!           expand potential Phi_n^j on a sphere-by-sphere basis [ needed for parallelism ]
!           ========================
            phiexp(ig,isph) = phi(ii)
!            
          endif
        enddo
      enddo
!
!     initialize forces
      fx = zero
!      
!     loop over atoms
      do isph = 1, nsph
!      
!       accumulate f += K_a contribution to < s , L' sigma >
        call fdoka( isph, sigma, xi(:,isph), basloc, dbsloc, vplm, vcos, vsin, fx(:,isph) ) 
!        
!       accumulate f += K_b contribution to < s , L' sigma >
        call fdokb( isph, sigma, xi,         basloc, dbsloc, vplm, vcos, vsin, fx(:,isph) ) 
!
!       accumulate f -= sum_n U_n^j' Phi_n^j xi(j,n) 
        call fdoga( isph,        xi,         phiexp,                           fx(:,isph) ) 
!        
      enddo

!!!      write(iout,1001)
!!! 1001 format(1x,'ddCOSMO forces [1] (atomic units):',/, &
!!!                1x,' atom',15x,'x',15x,'y',15x,'z')
!!!      do isph = 1, nsph
!!!        write(6,'(1x,i5,3f16.8)') isph, fx(:,isph)
!!!      enddo
!
!     deallocate workspaces
      deallocate( basloc, dbsloc, vplm, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'forces : [1] failed deallocation !'
      endif
!
!     allocate workspace
      allocate( zeta(ncav) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'forces : [3] failed allocation !'
      endif
!
!     initialize index
      ii = 0
!
!     loop over atoms
      do isph = 1, nsph
!
!       loop over gridpoints
        do ig = 1, ngrid
!
!         non-null contribution from grid point
          if ( ui(ig,isph).gt.zero ) then
!   
!           advance index
            ii = ii + 1
!
!           compute zeta(j,n)
!           =================
            zeta(ii) = xi(ig,isph)*w(ig)*ui(ig,isph)
!            
          endif
        enddo
      enddo
!
!     time computation of forces
      call system_clock( count = c2 )
!
!     printing
      if ( iprint.gt.0 ) then
!              
        write(iout,1010) dble(c2-c1)/dble(cr)
 1010   format(' computation time of the ddCOSMO forces = ',f8.3,' secs.')
! 
      endif
!
! --------------------------   modify here  --------------------------  
!
! here come the two contributions that depend on the solute, i.e., the electric field
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
! an example (efld). notice that the array csph contains the coordinates of the nuclei.
! the user will need to replace efld with his/her favorite routine.
!
!     allocate workspace
      allocate( ef(3,ncav) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'forces : [4] failed allocation !'
      endif
!
!     electric field produced by the charges, at the cavity points [ TARGETS ]
      call efld( n, charge, csph, ncav, ccav, ef )
!
!     initialize index
      ii = 0
!
!     loop over atoms
      do isph = 1, nsph
!
!       loop over gridpoints
        do ig = 1, ngrid
!
!         non-null contribution from integration point
          if ( ui(ig,isph).gt.zero ) then 
!
!           advance index
            ii = ii + 1
!
!           accumulate FIRST contribution to < s , g' >
            fx(:,isph) = fx(:,isph) + zeta(ii)*ef(:,ii)
!            
          endif
        enddo
      enddo
!
!!!      write(*,1002)
!!! 1002 format(1x,'ddCOSMO forces [2] (atomic units):',/, &
!!!                1x,' atom',15x,'x',15x,'y',15x,'z')
!!!      do isph = 1, nsph
!!!        write(6,'(1x,i5,3f16.8)') isph, fx(:,isph)
!!!      enddo

!     electric field produced by the cavity, at the nuclei [ TARGETS ]
      call efld( ncav, zeta, ccav, n, csph, ef )
!
!     loop over atoms
      do isph = 1, nsph
!      
!       accumulate SECOND contribution to < s , g' >
        fx(:,isph) = fx(:,isph) + ef(:,isph)*charge(isph)
!        
      enddo
      
!!!      write(*,1003)
!!! 1003 format(1x,'ddCOSMO forces [3] (atomic units):',/, &
!!!                1x,' atom',15x,'x',15x,'y',15x,'z')
!!!      do isph = 1, nsph
!!!        write(6,'(1x,i5,3f16.8)') isph, fx(:,isph)
!!!      enddo

!
! for point charges, there is no contribution from the derivatives of the psi vector.
! for quantum mechanical solutes, such a contribution needs to be handled via a numerical
! integration.
!
! --------------------------   end modify   --------------------------  
! 
!     deallocate workspaces
      deallocate( xi, phiexp, zeta, ef , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'forces : [2] failed deallocation !'
      endif
!
!     scale the forces time the cosmo factor
      fep = pt5*(eps-one)/eps
      fx  = fep*fx
!
!     printing      
      if ( iprint.ge.2 ) then
!              
        write(iout,1000)
 1000   format(1x,'ddCOSMO forces (atomic units):',/, &
                  1x,' atom',15x,'x',15x,'y',15x,'z')
!                  
        do isph = 1, nsph
          write(6,'(1x,i5,3f16.8)') isph, fx(:,isph)
        enddo
!        
      endif
!
!
endsubroutine forces
