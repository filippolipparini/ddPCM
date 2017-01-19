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
 1010   format(' the computation of the ddCOSMO part of the forces took ',f8.3,' seconds.')
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
!---------------------------------------------------------------------------------------
!
!
!
!
!---------------------------------------------------------------------------------------
subroutine check_derivativesCOSMO( )
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, zero, csph, rsph, memfree, &
                           ddinit, iquiet, basis, one, calcv2, fdoka, fdokb,     &
                           intrhs
!
      implicit none
      real*8 :: f(nbasis,nsph)
      real*8 :: e(nbasis,nsph)
      real*8 :: L(nbasis*nsph,nbasis*nsph),xi(ngrid,nbasis)
      real*8 :: dL(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: L_plus(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: xlm(nbasis),vplm(nbasis),vcos(lmax+1),vsin(lmax+1), &
                basloc(nbasis), pot(ngrid), dbsloc(3,nbasis)
      integer :: isph,jsph,i,j,ibeg,iend,nsph_save,icomp,ksph,iter
      real*8 :: eeps,err,rnorm
      integer, parameter :: niter = 6
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), s3(3)
      real*8 :: x(nsph), y(nsph), z(nsph), rwork(niter,nsph*3)
!
!---------------------------------------------------------------------------------------
!
!     construct L
!     -----------
      do isph = 1,nsph
!
        do jsph = 1,nsph
          do j = 1,nbasis
!
!           standard basis vector e_j
            e(:,   :) = zero
            e(j,jsph) = one
!
!           initialize
            pot = zero ; basloc = zero ; vplm = zero ; vcos = zero ; vsin = zero
!
!           compute L_i e_j
            call calcv2( .false., isph, pot, e, basloc, vplm, vcos, vsin )
!            
            ibeg = (isph-1)*nbasis+1
            iend = (isph-1)*nbasis+nbasis
            call intrhs( isph, pot, L( ibeg:iend , (jsph-1)*nbasis+j ) )
!
          enddo
        enddo
      enddo

!     construct dL
!     ------------
      do ksph = 1,nsph 
        do isph = 1,nsph
          do i = 1,nbasis
            do jsph = 1,nsph
              do j = 1,nbasis
!
!               standard basis vectors e_i, f_j
                e(:,:   )=zero
                f(:,:   )=zero
                e(i,isph)=one
                f(j,jsph)=one
!
!               auxiliary vector xi
                xi(:,:)=zero
                xi(:,isph)=basis(i,:)
!
!               initialize          
                s3 = zero ; dbsloc = zero
!
!               accumulate K_a contribution to < e , L' sigma >
                call fdoka( ksph, f, xi(:,isph), basloc, dbsloc, vplm, vcos, vsin, s3 ) 
!                
!               accumulate K_b contribution to < e , L' sigma >
                call fdokb( ksph, f, xi,         basloc, dbsloc, vplm, vcos, vsin, s3 ) 
!
!               store
                do icomp = 1,3
!
                  dL( (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp ) = s3(icomp)
!            
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
!     save initial DS
      nsph_save = nsph
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
!     set initial increment
      eeps=0.1d0
!
!     initialize
      rwork=zero
!
!     loop over increments
      do iter = 1,niter
!        
!       compute L^+
        do ksph = 1,nsph_save
!
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph) + eeps
            case(2) ; y(ksph) = y_save(ksph) + eeps
            case(3) ; z(ksph) = z_save(ksph) + eeps
            endselect
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!
!           initialize
            err=zero ; rnorm=zero
            L_plus(:,:,ksph,icomp)=zero
!
!           build L^+
            do isph = 1,nsph_save
!            
              do jsph = 1,nsph_save
                do j = 1,nbasis
!
!                 standard basis vector e_j
                  e(:,:   )=zero
                  e(j,jsph)=one
!                  
!                 initialize
                  pot = zero ; basloc = zero ; vplm = zero ; vcos = zero ; vsin = zero
!
!                 compute L^+_i e_j
                  call calcv2( .false., isph, pot, e, basloc, vplm, vcos, vsin )
!
                  ibeg = (isph-1)*nbasis+1
                  iend = (isph-1)*nbasis+nbasis
                  call intrhs( isph, pot, L_plus( ibeg:iend , (jsph-1)*nbasis+j , ksph , icomp ) )
!
!                 accumulate error
                  do i = 1,nbasis
!                  
                    err   = err   + ( ( L_plus((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) -           &
                                        L(     (isph-1)*nbasis+i,(jsph-1)*nbasis+j           ) ) / eeps -  &
                                        dL(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)                 )**2
!                                        
                    rnorm = rnorm + (   dL(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)                 )**2
!
                  enddo

                enddo
              enddo
            enddo
!            
!!!!           numerical derivatives
!!!            write(*,1007) ksph,icomp
!!! 1007       format(' A(r+r_',i2,',',i1,') =')
!!!! 
!!!            do isph = 1,nsph
!!!              do i = 1,nbasis
!!!                write(*,"(4x,300(e12.5,2x))") ( A_plus((isph-1)*nbasis+i,j,ksph,icomp) , j=1,nbasis*nsph )
!!!              enddo
!!!            enddo
!!!            write(*,*)''
!!!!            
!!!            write(*,1008) ksph,icomp
!!! 1008       format('( A(r+r_',i2,',',i1,') - A(r) ) / eps =')
!!!! 
!!!            do isph = 1,nsph
!!!              do i = 1,nbasis
!!!                write(*,"(4x,300(e12.5,2x))") ( ( A_plus((isph-1)*nbasis+i,j,ksph,icomp) - &
!!!                                                     A((isph-1)*nbasis+i,j))/eeps , j=1,nbasis*nsph )
!!!              enddo
!!!            enddo
!!!            write(*,*)''
!
!           store error
            err=sqrt(err) ; rnorm=sqrt(rnorm)
!
            if ( rnorm.gt.1.E-12 ) then
!
              rwork(iter,(ksph-1)*3+icomp) = sqrt( err / rnorm )
!              
            endif
!
          enddo
        enddo
!
!       update increment
        eeps = eeps/2.d0
!        
      enddo

      eeps = eeps*2.d0
!
!     printing
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dL / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!!!!     printing
!!!      do ksph = 1,nsph
!!!        do icomp = 1,3
!!!!        
!!!!         analytical derivatives
!!!          write(*,1005) ksph,icomp
!!! 1005     format(' dA / dr_',i2,',',i1,' =')
!!!!
!!!          do isph = 1,nsph
!!!            do i = 1,nbasis
!!!              write(*,"(4x,300(e12.5,2x))") ( dA((isph-1)*nbasis+i,j,ksph,icomp), j=1,nbasis*nsph )
!!!            enddo
!!!          enddo
!!!!          
!!!!         numerical derivatives
!!!          write(*,1006) ksph,icomp
!!! 1006     format(' ( A(r+r_',i2,',',i1,') - A(r) ) / eps =')
!!!! 
!!!          do isph = 1,nsph
!!!            do i = 1,nbasis
!!!              write(*,"(4x,300(e12.5,2x))") &
!!!              ( ( A_plus((isph-1)*nbasis+i,j,ksph,icomp)- &
!!!                       A((isph-1)*nbasis+i,j)             ) / eeps , j=1,nbasis*nsph )
!!!            enddo
!!!          enddo
!!!          write(*,*)''
!!!! 
!!!        enddo
!!!      enddo
!
!     restore DS
      call memfree
      call ddinit( nsph_save, x_save, y_save, z_save, r_save )
!
!     deactivate quiet flag      
      iquiet = .false.
!
!
endsubroutine check_derivativesCOSMO
!---------------------------------------------------------------------------------------
!
!
!
!
!---------------------------------------------------------------------------------------
subroutine check_forcesCOSMO( E0, charge, f )
!
      use ddcosmo , only : nbasis, nsph, iquiet, csph, rsph, memfree, ddinit, &
                           ncav, ccav, ngrid, zero, itsolv
!                           
      implicit none
      real*8,                    intent(in) :: E0
      real*8, dimension(  nsph), intent(in) :: charge
      real*8, dimension(3,nsph), intent(in) :: f
!
      integer,parameter :: niter = 6
!
      real*8 :: phi(ncav), psi(nbasis,nsph), sigma(nbasis,nsph)
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), phi_eps(nbasis,nsph)
      real*8 :: x(nsph), y(nsph), z(nsph), rwork(niter,nsph*3),rrate(niter,nsph*3)
      real*8 :: E_plus, err, eeps
      integer :: iter, icomp, ksph, nsph_save, j
!
!---------------------------------------------------------------------------------------
!
!     activate quiet flag
      iquiet = .true.
!
!     save initial DS
      nsph_save = nsph
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
!     set initial increment
      eeps=0.1d0
!
!     initialize
      rwork = zero ; rrate = zero
!
!     loop over increments
      do iter = 1,niter
!        
!       loop over d / dr_k
        do ksph = 1,nsph_save
!
!         loop over components of d / dr_k
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph) + eeps
            case(2) ; y(ksph) = y_save(ksph) + eeps
            case(3) ; z(ksph) = z_save(ksph) + eeps
            endselect
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!            
!           potential Phi and Psi vector
            call mkrhs( nsph, charge, x, y, z, ncav, ccav, phi, nbasis, psi )
!
!           solve COSMO equations       
            E_plus = zero ; sigma = zero
            call itsolv( .false., phi, psi, sigma, E_plus )
!
            if ( abs( f(icomp,ksph) ).gt.1.E-12 ) then
!                    
!             compute relative error
              err = abs( (E_plus - E0) / eeps + f(icomp,ksph) ) / abs( f(icomp,ksph) )
!
!             store
              rwork(iter,(ksph-1)*3+icomp) = err
!              
!             compute rate
              if ( iter.gt.1 ) then 
                rrate(iter,(ksph-1)*3+icomp) =  log( rwork(iter-1,(ksph-1)*3+icomp) / &
                                                     rwork(iter  ,(ksph-1)*3+icomp)   ) / log(0.5d0)  
              endif
            endif
!
          enddo
        enddo
!
        eeps = eeps / 2.d0
!
      enddo
!      
!     print relative error
      write(*,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!     print rate of convergence
      write(*,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!      
!     restore DS
      call memfree
      call ddinit( nsph_save, x_save, y_save, z_save, r_save )
!
!     deactivate quiet flag      
      iquiet = .false.
!
!
endsubroutine check_forcesCOSMO
