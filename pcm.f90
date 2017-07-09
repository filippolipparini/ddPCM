subroutine pcm(star, cart, phi, glm, phi_eps)
use ddcosmo
implicit none
!
! wrapper for the linear solvers for IEFPCM. The IEFPCM equation we are want to solve is
!
!     R_\eps \Phi_\eps = R_\infty \Phi.
! 
! the right-hand side is therefore g = R_\infty \Phi.
!
! input:
! 
!   star     logical, true:  solve the adjoint PCM equations,
!                     false: solve the PCM equatiosn
!
!   cart     logical, true:  the right-hand side for the PCM has to be assembled 
!                            inside this routine and the unscaled potential at the 
!                            external points of the cavity is provided in phi. 
!                     false: the right-hand side for the PCM equations is provided
!                            in glm.
!                     cart is not referenced if star is true. 
!
!   phi      real,    contains the potential at the external cavity points if star is
!                     false and cart is true.
!                     phi is not referenced in any other case.
!
!   glm      real,    contains the right-hand side for the PCM (adjoint) equations. 
!                     if star is false and cart is true, glm is not referenced. 
!
! output:
!
!   phi_eps: real,    the solution to the COSMO (adjoint) equations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! this routine performs the following operations:
!
!   - allocates memory for the linear solvers, and fixes dodiag.
!     This parameters controls whether the diagonal part of the matrix is considered 
!     in matvec, which depends on the solver used. It is false for jacobi_diis and
!     true for GMRES. 
!
!   - if star is false and cart is true, assembles the right-hand side for the PCM
!     equations. Note that for GMRES, a preconditioner is applied.
!
!   - computes a guess for the solution (using the preconditioner)
!
!   - calls the required iterative solver
!
logical,                         intent(in)    :: star, cart
real*8,  dimension(ncav),        intent(in)    :: phi
real*8,  dimension(nbasis,nsph), intent(in)    :: glm
real*8,  dimension(nbasis,nsph), intent(inout) :: phi_eps
!
integer              :: isph, istatus, n_iter, info, c1, c2, cr
real*8               :: tol, r_norm
logical              :: ok
!
real*8, allocatable  :: g(:,:), rhs(:,:), work(:,:), x(:,:), u(:), ulm(:), basloc(:), &
                        vplm(:), vcos(:), vsin(:)
!
integer, parameter   :: gmm = 20, gmj = 25
!
external             :: rx, prx, precx, hnorm, rstarx, prstarx
!
! set a few parameters for the solver and matvec routine:
!
tol     = 10.0d0**(-iconv)
n_iter  = 300
!
! initialize the timer:
!
call system_clock(count_rate=cr)
call system_clock(count=c1)
!
! set solver-specific options:
!
if (isolver .eq. 0) then
!
! jacobi_diis
!
  do_diag = .false.
!
else
!
! GMRES
!
  do_diag = .true.
!
  allocate (work(nsph*nbasis,0:2*gmj+gmm+2 -1), stat=istatus)
  if (istatus .ne. 0) then
    write(*,*) ' pcm: [1] failed allocation for GMRES'
    stop
  end if
  work  = zero
!
end if
!
if (.not. star) then
!
! solve LX = g.
!
  allocate (rhs(nbasis,nsph), stat=istatus)
!
  rhs = zero
!
  if (istatus .ne. 0) then
    write(*,*) ' pcm: [2] failed allocation'
  end if
!
  if (cart) then
    allocate (g(ngrid,nsph), x(nbasis,nsph), u(ngrid), ulm(nbasis), basloc(nbasis), &
              vplm(nbasis), vcos(lmax+1), vsin(lmax+1), stat=istatus)
!
    if (istatus .ne. 0) then
      write(*,*) ' pcm: [3] failed allocation'
    end if
!
!   we need to assemble the right-hand side by weighting the potential
!   and calling intrhs. Start weighting the potential...
! 
    call wghpot(phi, g)
!
!   and compute its multipolar expansion
!
    do isph = 1, nsph
      call intrhs(isph, g(:,isph), x(:,isph))
    end do
!
!   now, apply R_\infty:
!
    do isph = 1, nsph
      call mkrvec(isph, zero, x, rhs(:,isph), ulm, u, basloc, vplm, vcos, vsin)
    end do
!
    if (isolver .eq. 1) then
      do_diag = .true.
      x = rhs
      call precx(nbasis*nsph, x, rhs)
    else
      do_diag = .false.
    end if
    deallocate (g, x, u, ulm, basloc, vplm, vcos, vsin)
!
  else
!
    if (isolver .eq. 0) then
      rhs = glm
    else if (isolver .eq. 1) then
!     rhs = glm
      call precx(nbasis*nsph, glm, rhs)
    end if
  end if
!
! assemble a guess:
!
  call precx(nsph*nbasis,rhs,phi_eps)
!
! call the solver:
!
  if (isolver .eq. 0) then
!
!   jacobi/diis
!
    call jacobi_diis(nsph*nbasis, iprint, ndiis, 3, tol, rhs, phi_eps, n_iter, ok, rx, precx)
    do_diag = .true.
!
  else if (isolver .eq. 1) then
!
!   gmres. the gmres solver can not handle preconditioners, so we will solve 
!  
!     PLX = Pg,
!
!   where P is a jacobi preconditioner. note thus the plx matrix-vector multiplication routine.
!
    call gmresr(iprint.gt.0, nsph*nbasis, gmj, gmm, rhs, phi_eps, work, tol, 'abs', n_iter, r_norm, prx, info)
    ok = info .eq. 0
!
  end if
!
! esolv = pt5 * ((eps - one)/eps) * sprod(nsph*nbasis,sigma,psi)
!
  deallocate (rhs)
!
else
!
  allocate (rhs(nbasis,nsph), stat=istatus)
  if (istatus .ne. 0) then
    write(6,*) 'pcm: [4] allocation failed '
    stop
  end if
!
  if (isolver .eq. 0) then
    rhs = glm
  else if (isolver .eq. 1) then
    call precx(nbasis*nsph,glm,rhs)
  end if
!
! assemble a guess:
!
  call precx(nbasis*nsph,rhs,phi_eps)
!
! call the solver
!
  if (isolver .eq. 0) then
!
!   jacobi/diis
!
    call jacobi_diis(nsph*nbasis, iprint, ndiis, 3, tol, rhs, phi_eps, n_iter, ok, rstarx, precx)
!fl
    do_diag = .true.
!
  else if (isolver .eq. 1) then
!
!   gmres. the gmres solver can not handle preconditioners, so we will solve 
!  
!     PR_\eps^* Phi_eps = P g,
!
!   where P is a jacobi preconditioner. note thus the plx matrix-vector multiplication routine.
!
    call gmresr(iprint.gt.0, nsph*nbasis, gmj, gmm, rhs, phi_eps, work, tol, 'abs', n_iter, r_norm, prstarx, info)
    ok = info .eq. 0
!
  end if
!fl   do isph = 1, nsph
!fl     sigma(:,isph) = facl(:)*psi(:,isph)
!fl   end do
!fl !
!fl ! call the solver:
!fl !
!fl   if (isolver .eq. 0) then
!fl     call jacobi_diis(nsph*nbasis, iprint, ndiis, 4, tol, psi, sigma, n_iter, ok, lstarx, ldm1x, hnorm)
!fl   else if (isolver .eq. 1) then
!fl     allocate (rhs(nbasis,nsph), stat=istatus)
!fl     if (istatus .ne. 0) then
!fl       write(*,*) 'cosmo: [3] failed allocation'
!fl       stop
!fl     end if
!fl !
!fl !   gmres. the gmres solver can not handle preconditioners, so we will solve 
!fl !  
!fl !     PL*S = P\Psi,
!fl !
!fl !   where P is a jacobi preconditioner. note thus the pstarlx matrix-vector multiplication routine.
!fl !
!fl     call ldm1x(nsph*nbasis,psi,rhs)
!fl     call gmresr(iprint.gt.0, nsph*nbasis, gmj, gmm, rhs, sigma, work, tol, 'abs', n_iter, r_norm, plstarx, info)
!fl     ok = info .eq. 0
!fl     deallocate (rhs)
!fl   end if
!
end if
!
if (isolver .eq. 1) deallocate (work)
!
if (.not. ok) then
  if (star) then
    write(iout,1020) 'adjoint '
  else
    write(iout,1020) ''
  end if
  stop
end if
!
call system_clock(count=c2)
!
if (iprint.gt.0) then
  if (star) then
    write(iout,1010) 'adjoint ', dble(c2-c1)/dble(cr)
  else
    write(iout,1010) '', dble(c2-c1)/dble(cr)
  end if
  write(iout,*)
end if
!
 1010 format(' the solution to the ddCOSMO ',a,'equations took ',f8.3,' seconds.')
 1020 format(' ddCOSMO did not converge! Aborting...')
!
return
end
