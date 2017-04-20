subroutine numgrad(fx,x,y,z,rvdw,charge)
use ddcosmo
implicit none
real*8,  intent(inout), dimension(nsph) :: x, y, z, rvdw, charge
real*8,  intent(in), dimension(3,nsph)  :: fx
!
integer :: i, k, n
real*8  :: eph, emh, h
!
real*8, allocatable :: phi(:), psi(:,:), g(:,:), sigma(:,:), fnum(:,:)
!
allocate (phi(2*ncav), psi(nbasis,nsph), g(ngrid,nsph), sigma(nbasis,nsph), fnum(3,nsph))
!
write(6,*) 'h =?'
read(5,*) h
!
do i = 1, nsph
  x(i) = x(i) + h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call itsolv(.false.,phi,psi,sigma,eph)
  call memfree
  x(i) = x(i) - two*h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call itsolv(.false.,phi,psi,sigma,emh)
  call memfree
  x(i) = x(i) + h
  fnum(1,i) = (eph - emh)/(two*h)
  y(i) = y(i) + h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call itsolv(.false.,phi,psi,sigma,eph)
  call memfree
  y(i) = y(i) - two*h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call itsolv(.false.,phi,psi,sigma,emh)
  call memfree
  y(i) = y(i) + h
  fnum(2,i) = (eph - emh)/(two*h)
  z(i) = z(i) + h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call itsolv(.false.,phi,psi,sigma,eph)
  call memfree
  z(i) = z(i) - two*h
  call ddinit(nsph,x,y,z,rvdw)
  call mkrhs(nsph,charge,x,y,z,ncav,ccav,phi,nbasis,psi)
  call itsolv(.false.,phi,psi,sigma,emh)
  call memfree
  z(i) = z(i) + h
  fnum(3,i) = (eph - emh)/(two*h)
end do
!
do i = 1, nsph
  write(6,'(i4,6f16.8)') i, (fx(k,i),fnum(k,i), k = 1, 3)
end do
!
stop
return
end
