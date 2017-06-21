subroutine cosmo(star, phi, psi, sigma, esolv)
use ddcosmo
implicit none
logical,                         intent(in)    :: star
real*8,  dimension(ncav),        intent(in)    :: phi
real*8,  dimension(nbasis,nsph), intent(in)    :: psi
real*8,  dimension(nbasis,nsph), intent(inout) :: sigma
real*8,                          intent(inout) :: esolv
!
integer :: isph, istatus, n_iter
real*8  :: tol
logical :: ok
!
external :: lx, ldm1x, hnorm, lstarx
!
real*8, allocatable :: g(:,:), rhs(:,:)
!
if (.not. star) then
  allocate (g(ngrid,nsph), rhs(nbasis,nsph), stat=istatus)
  if (istatus .ne. 0) then
    write(*,*) ' cosmo: [1] failed allocation'
  end if
!
! weight the potential
!
  call wghpot(phi, g)
!
! and compute its multipolar expansion
!
  do isph = 1, nsph
    call intrhs(isph, g(:,isph), rhs(:,isph))
  end do
!
! set a few parameters for the solver and matvec routine:
!
  do_diag = .false.
  tol     = 10.0d0**(-iconv)
  n_iter  = 100
!
! assemble a guess:
!
  do isph = 1, nsph
    sigma(:,isph) = facl(:)*rhs(:,isph)
  end do
!
! call the solver:
!
  call jacobi_diis(nsph*nbasis, iprint, ndiis, 4, tol, rhs, sigma, n_iter, ok, lx, ldm1x, hnorm)
!
  esolv = pt5 * ((eps - one)/eps) * sprod(nsph*nbasis,sigma,psi)
!
  deallocate (g, rhs)
!
else
!
! set a few parameters for the solver and matvec routine:
!
  do_diag = .false.
  tol     = 10.0d0**(-iconv) 
  n_iter  = 100
!
! assemble a guess:
!
  do isph = 1, nsph
    sigma(:,isph) = facl(:)*psi(:,isph)
  end do
!
! call the solver:
!
  call jacobi_diis(nsph*nbasis, iprint, ndiis, 4, tol, psi, sigma, n_iter, ok, lstarx, ldm1x, hnorm)
!
end if
!
return
end
