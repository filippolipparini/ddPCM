subroutine cosmo( star, cart, phi, glm, psi, sigma, esolv )
!
!---------------------------------------------------------------------------------------
! Purpose : wrapper for the linear solvers for COSMO equation
!
!             L sigma = G 
!
!           and adjoint COSMO equation
!
!             L^* sigma = Psi
!
!---------------------------------------------------------------------------------------
!
! input:
! 
!   star   logical, true:  solve the adjoint COSMO equations,
!                   false: solve the COSMO equatiosn
!
!   cart   logical, true:  the right-hand side for the COSMO has to be assembled 
!                          inside this routine and the unscaled potential at the 
!                          external points of the cavity is provided in phi. 
!                   false: the right-hand side for the COSMO equations is provided
!                          in glm.
!                   cart is not referenced if star is true. 
!
!   phi    real,    contains the potential at the external cavity points if star is
!                   false and cart is true.
!                   phi is not referenced in any other case.
!
!   glm    real,    contains the right-hand side for the COSMO equations if star is
!                   false and cart is false.
!                   glm is not referenced in any other case
!
!   psi    real,    the psi vector. it is used to compute the energy if star is false,
!                   as a right-hand side if star is true.
!
! output:
!
!   sigma: real,    the solution to the COSMO (adjoint) equations
!
!   esolv: real,    if star is false, the solvation energy.
!                   if star is true, it is not referenced.
!            
!---------------------------------------------------------------------------------------
! This routine performs the following operations :
!
!   - allocates memory for the linear solvers
!
!   - if star is false and cart is true, assembles the right-hand side for the COSMO
!     equations. 
!
!   - computes a guess for the solution (using the inverse diagonal);
!
!   - calls the iterative solver;
!
!   - if star is false, computes the solvation energy.
!---------------------------------------------------------------------------------------
!
      use ddcosmo , only : iprint, ncav, nylm, nsph, iconv, zero, ngrid, ndiis,      &
                           wghpot, intrhs, facl, pt5, eps, sprod, iout, one, prtsph
!      
      implicit none
      logical,                         intent(in)    :: star, cart
      real*8,  dimension(ncav),        intent(in)    :: phi
      real*8,  dimension(nylm,nsph), intent(in)      :: glm, psi
      real*8,  dimension(nylm,nsph), intent(inout)   :: sigma
      real*8,                          intent(inout) :: esolv
!
      integer              :: isph, istatus, n_iter, info, c1, c2, cr
      real*8               :: tol, r_norm
      logical              :: ok
!
      real*8, allocatable  :: g(:,:), rhs(:,:), work(:,:)
!
      external             :: lx, ldm1x, hnorm, lstarx
!
!---------------------------------------------------------------------------------------
!
!     parameters for the solver and matvec routine
      tol     = 10.0d0**(-iconv)
      n_iter  = 200
!
!     initialize the timer
      call system_clock(count_rate=cr)
      call system_clock(count=c1)
!
!     DIRECT COSMO EQUATION L X = g
!     =============================
!
      if ( .not.star ) then
!
!       allocate workspace for rhs
        allocate( rhs(nylm,nsph) , stat=istatus )
        if (istatus .ne. 0) then
          write(*,*) ' cosmo: [2] failed allocation'
        endif
!
!       1. RHS
!       ------

!       assemble rhs
        if ( cart ) then
!
!         allocate workspace for weighted potential
          allocate( g(ngrid,nsph) , stat=istatus )
          if (istatus .ne. 0) then
            write(*,*) ' cosmo: [3] failed allocation'
          endif
!
!         weight the potential...
          call wghpot( phi, g )
!
!         ... and compute its multipolar expansion
          do isph = 1, nsph
            call intrhs( isph, g(:,isph), rhs(:,isph) )
          enddo
!
!         deallocate workspace
          deallocate( g , stat=istatus )
          if ( istatus.ne.0 ) then
            write(*,*) 'cosmo: [1] failed deallocation'
          endif
!
!       no need to manipulate rhs
        else
!
          rhs = glm
!
        endif
!
        if (iprint.ge.4) call prtsph('rhs of the ddCOSMO equation',nsph,0,rhs)
!
!       2. INITIAL GUESS
!       ----------------
!
        do isph = 1, nsph
          sigma(:,isph) = facl(:)*rhs(:,isph)
        enddo
!
!       3. SOLVER CALL
!       --------------
!
!       Jacobi method : 
!
!         L X = ( diag + offdiag ) X = g   ==>    X = diag^-1 ( g - offdiag X_guess )
!
!         action of  diag^-1 :  ldm1x
!         action of  offdiag :  lx
!
          call jacobi_diis( nsph*nylm, iprint, ndiis, 4, tol, rhs, sigma, n_iter, ok, lx, ldm1x, hnorm )
!
!       4. SOLVATION ENERGY
!       -------------------
!
        esolv = pt5 * ((eps - one)/eps) * sprod( nsph*nylm, sigma, psi )
!
!       deallocate workspace
        deallocate( rhs , stat=istatus )
        if ( istatus.ne.0 ) then
           write(*,*) 'cosmo: [2] failed deallocation'
        endif
!
!
!     ADJOINT COSMO EQUATION L^* X = Psi
!     ==================================
!
      else
!
        if (iprint.ge.4) call prtsph('rhs of the ddCOSMO adjoint equation',nsph,0,psi)
!
!       1. INITIAL GUESS
!       ----------------
!
        do isph = 1, nsph
          sigma(:,isph) = facl(:)*psi(:,isph)
        enddo
!
!       2. SOLVER CALL
!       --------------
!
!       Jacobi method : see above
!
        call jacobi_diis( nsph*nylm, iprint, ndiis, 4, tol, psi, sigma, n_iter, ok, lstarx, ldm1x, hnorm )
!          
      endif
!
!     check solution
      if ( .not.ok ) then
!              
        if ( star ) then
          write(iout,1020)
 1020     format(' adjoint ddCOSMO did not converge! Aborting...')
        else
          write(iout,1021)
 1021     format(' ddCOSMO did not converge! Aborting...')
        endif
!        
        stop
!        
      endif
!
      call system_clock(count=c2)
!
!     printing
      if ( iprint.gt.0 ) then
!              
        write(iout,*)
!        
!       adjoint
        if ( star ) then
          write(iout,1010) dble(c2-c1)/dble(cr)
 1010     format(' solution time of ddCOSMO adjoint equations L^* \sigma = ',f8.3,' secs.')
!
!       direct
        else
          write(iout,1011) dble(c2-c1)/dble(cr)
 1011     format(' solution time of ddCOSMO direct equations L \sigma = ',f8.3,' secs.')
! 
        endif
!        
        write(iout,*)
! 
      endif
!
!
endsubroutine cosmo
