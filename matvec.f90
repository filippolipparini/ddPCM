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
