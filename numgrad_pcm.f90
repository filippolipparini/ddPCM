subroutine numgrad_pcm(fx,x,y,z,rvdw,charge)
use ddcosmo
implicit none
real*8,  intent(inout), dimension(nsph) :: x, y, z, rvdw, charge
real*8,  intent(in), dimension(nsph,3)  :: fx
!
integer :: i, k, n
real*8  :: eph, emh, h, xx(1)
!
real*8, allocatable :: phi(:), psi(:,:), g(:,:), sigma(:,:), fnum(:,:), phi_eps(:,:)
!
allocate (phi(2*ncav), psi(nbasis,nsph), g(ngrid,nsph), sigma(nbasis,nsph), &
          phi_eps(nbasis,nsph), fnum(3,nsph))
!
write(6,*) 'h =?'
read(5,*) h
!
do i = 1, nsph
  x(i) = x(i) + h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call pcm(.false., .true., .true., phi, xx, phi_eps)
  call cosmo(.false., .false., xx, phi_eps, psi, sigma, eph)
  call memfree
  x(i) = x(i) - two*h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call pcm(.false., .true., .true., phi, xx, phi_eps)
  call cosmo(.false., .false., xx, phi_eps, psi, sigma, emh)
  call memfree
  x(i) = x(i) + h
  fnum(1,i) = (eph - emh)/(two*h)
  y(i) = y(i) + h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call pcm(.false., .true., .true., phi, xx, phi_eps)
  call cosmo(.false., .false., xx, phi_eps, psi, sigma, eph)
  call memfree
  y(i) = y(i) - two*h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call pcm(.false., .true., .true., phi, xx, phi_eps)
  call cosmo(.false., .false., xx, phi_eps, psi, sigma, emh)
  call memfree
  y(i) = y(i) + h
  fnum(2,i) = (eph - emh)/(two*h)
  z(i) = z(i) + h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call pcm(.false., .true., .true., phi, xx, phi_eps)
  call cosmo(.false., .false., xx, phi_eps, psi, sigma, eph)
  call memfree
  z(i) = z(i) - two*h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call pcm(.false., .true., .true., phi, xx, phi_eps)
  call cosmo(.false., .false., xx, phi_eps, psi, sigma, emh)
  call memfree
  z(i) = z(i) + h
  fnum(3,i) = (eph - emh)/(two*h)
end do
!
do i = 1, nsph
  write(6,'(i4,6f16.8)') i, (fx(i,k),fnum(k,i), k = 1, 3)
end do
!
stop
return
end
