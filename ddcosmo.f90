module ddcosmo
implicit none
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
integer, parameter :: ndiis=25, iout=6, nngmax=100
real*8,  parameter :: zero=0.d0, pt5=0.5d0, one=1.d0, two=2.d0, four=4.d0
!
integer :: nsph, ngrid, ncav, lmax, nbasis, iconv, igrad, &
           iprint, nproc, memuse, memmax
real*8  :: eps, eta, pi, sq2
logical :: grad
!
integer, allocatable :: inl(:), nl(:)
real*8,  allocatable :: rsph(:), csph(:,:), ccav(:,:)
real*8,  allocatable :: w(:), grid(:,:), basis(:,:)
real*8,  allocatable :: fact(:), facl(:), facs(:)
real*8,  allocatable :: fi(:,:), ui(:,:), zi(:,:,:)
!
contains
subroutine ddinit(n,x,y,z,rvdw)
implicit none
!
! allocate the various arrays needed for ddcosmo,
! assemble the cavity and the various associated geometrical quantities.
!
integer,               intent(in) :: n
real*8,  dimension(n), intent(in) :: x, y, z, rvdw
!
integer :: isph, jsph, i, ii, lnl, l, ind, m, igrid, inear, jnear
real*8  :: fac, fl, ffl, fnorm, d2, r2, v(3), vv, t, xt, swthr
!
real*8,  allocatable :: vcos(:), vsin(:), vplm(:)
integer, parameter   :: nllg=32
!
integer, dimension(nllg) :: ng0
data ng0/6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974, &
         1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810/
!
! openMP parallelization:
!
if (nproc.eq.0) nproc = 1
!$ call omp_set_num_threads(nproc)

pi   = four*atan(one)
sq2  = sqrt(two)
nsph = n
!
! choose the lebedev grid with number of points closest to ngrid:
!
igrid = 0
inear = 100000
do i = 1, nllg
  jnear = iabs(ng0(i)-ngrid)
  if (jnear.lt.inear) then
    inear = jnear
    igrid = i
  end if
end do
ngrid = ng0(igrid)
!
! print a nice header:
!
call header
!
! allocate:
!
grad   = igrad.ne.0
nbasis = (lmax+1)*(lmax+1)
allocate (rsph(nsph),csph(3,nsph))
allocate (w(ngrid),grid(3,ngrid),basis(nbasis,ngrid))
allocate (inl(nsph+1),nl(nsph*nngmax))
allocate (fi(ngrid,nsph),ui(ngrid,nsph))
if (grad) allocate(zi(3,ngrid,nsph))
allocate (fact(2*lmax+1),facl(nbasis),facs(nbasis))
!
memuse = memuse + 4*nsph + 4*ngrid + nbasis*ngrid + nsph+1 + nsph*nngmax + &
         2*ngrid*nsph + 2*lmax+1 + 2*nbasis
if (grad) memuse = memuse + 3*ngrid*nsph
memmax = max(memmax,memuse)
!
! precompute various quantities (diagonal blocks of ddcosmo matrix,
! normalization factors for spherical harmonics, factorials)
!
fact(1) = one
fact(2) = one
do i = 3, 2*lmax + 1
  fact(i) = dble(i-1)*fact(i-1)
end do
!
do l = 0, lmax
  ind = l*l + l + 1
  fl  = (two*dble(l) + one)/(four*pi)
  ffl = sqrt(fl)
  facl(ind-l:ind+l) = fl
  facs(ind) = ffl
  do m = 1, l
    fnorm = sq2*ffl*sqrt(fact(l-m+1)/fact(l+m+1))
    if (mod(m,2).eq.1) fnorm = -fnorm
    facs(ind+m) = fnorm
    facs(ind-m) = fnorm
  end do
end do
!
! set the centers and radii of the spheres:
!
csph(1,:) = x
csph(2,:) = y
csph(3,:) = z
rsph      = rvdw
!
! load a lebedev grid:
!
call llgrid(ngrid,w,grid)
!
! build a basis of spherical harmonics at the gridpoints:
!
allocate (vplm(nbasis),vcos(lmax+1),vsin(lmax+1))
memuse = memuse + nproc*(nbasis + 2*lmax + 2)
memmax = max(memmax,memuse)
!$omp parallel do default(shared) private(i,vplm,vcos,vsin)
do i = 1, ngrid
  call ylmbas(grid(:,i),basis(:,i),vplm,vcos,vsin)
end do
!$omp end parallel do
deallocate (vplm,vcos,vsin)
!
memuse = memuse - nproc*(nbasis + 2*lmax + 2)
!
! build a neighbor list:
!
ii  = 1
lnl = 0
do isph = 1, nsph
  inl(isph) = lnl + 1
  do jsph = 1, nsph
    if (isph.ne.jsph) then
      d2 = (csph(1,isph) - csph(1,jsph))**2 + (csph(2,isph) - csph(2,jsph))**2 + (csph(3,isph) - csph(3,jsph))**2
      r2 = (rsph(isph) + rsph(jsph))**2
      if (d2.le.r2) then
        nl(ii) = jsph
        ii  = ii + 1
        lnl = lnl + 1
      end if
    end if
  end do
end do
inl(nsph+1) = lnl+1
!
! build ui, fi and zi:
!
fi = zero
ui = zero
if (grad) zi = zero
!$omp parallel do default(shared) private(isph,i,ii,jsph,v,vv,t,xt,swthr,fac)
do isph = 1, nsph
  do i = 1, ngrid
    do ii = inl(isph), inl(isph+1) - 1
      jsph = nl(ii)
      v(:)  = csph(:,isph) + rsph(isph)*grid(:,i) - csph(:,jsph)
      vv = sqrt(dot_product(v,v))
      t  = vv/rsph(jsph)
      xt  = fsw(t,eta*rsph(jsph))
      swthr  = one - eta*rsph(jsph)
      if (grad .and. (t.lt.one .and. t.gt.swthr)) then
        fac = dfsw(t,eta*rsph(jsph))/rsph(jsph)
        zi(:,i,isph) = zi(:,i,isph) + fac*v(:)/vv
      end if
      fi(i,isph) = fi(i,isph) + xt
    end do
    if (fi(i,isph).le.one) ui(i,isph) = one - fi(i,isph)
  end do
end do
!$omp end parallel do
!
! count the number of external points and allocate the cavity:
!
ncav = 0
do isph = 1, nsph
  do i = 1, ngrid
    if (ui(i,isph).gt.zero) ncav = ncav + 1
  end do
end do
!
allocate (ccav(3,ncav))
memuse = memuse + 3*ncav
memmax = max(memmax,memuse)
!
! fill ccav with the coordinates of the external points:
!
ii = 0
do isph = 1, nsph
  do i = 1, ngrid
    if (ui(i,isph).gt.zero) then
      ii = ii + 1
      ccav(:,ii) = csph(:,isph) + rsph(isph)*grid(:,i)
    end if
  end do
end do
return
end subroutine ddinit
!
subroutine memfree
implicit none
!
! deallocate the arrays:
!
if(allocated(rsph))  deallocate(rsph)  
if(allocated(csph))  deallocate(csph)  
if(allocated(ccav))  deallocate(ccav)  
if(allocated(w))     deallocate(w)     
if(allocated(grid))  deallocate(grid)  
if(allocated(basis)) deallocate(basis) 
if(allocated(inl))   deallocate(inl)   
if(allocated(nl))    deallocate(nl)    
if(allocated(fact))  deallocate(fact)  
if(allocated(facl))  deallocate(facl)  
if(allocated(facs))  deallocate(facs)  
if(allocated(ui))    deallocate(ui)
if(allocated(fi))    deallocate(fi)
if(allocated(zi))    deallocate(zi)
!
memuse = memuse - 4*nsph - 4*ngrid - nbasis*ngrid - nsph-1 - nsph*nngmax - &
         2*ngrid*nsph - 2*lmax-1 - 2*nbasis
if (grad) memuse = memuse - 3*ngrid*nsph
end subroutine memfree
!
real*8 function sprod(n,u,v)
implicit none
integer,               intent(in) :: n
real*8,  dimension(n), intent(in) :: u, v
!
integer :: i
real*8  :: ss
!
ss = zero
do i = 1, n
  ss = ss + u(i)*v(i)
end do
sprod = ss
return
end function sprod
!
real*8 function fsw(t,eta)
implicit none
real*8, intent(in) :: t, eta
!
real*8 :: a, b, flow
real*8, parameter :: f6=6.0d0, f10=10.d0, f12=12.d0, f15=15.d0
flow = one - eta
if (t.ge.one) then
  fsw = zero
else if (t.le.flow) then
  fsw = one
else
  a = f15*eta - f12
  b = f10*eta*eta - f15*eta + f6
  fsw = ((t-one)*(t-one)*(one-t)*(f6*t*t + a*t + b))/(eta**5)
end if
return
end function fsw
!
real*8 function dfsw(t,eta)
implicit none
real*8, intent(in) :: t, eta
!
! switching function derivative for ddCOSMO regularization.
!
real*8  flow
real*8, parameter :: f30=30.0d0
!
flow = one - eta
if (t.ge.one) then
  dfsw = zero
else if (t.le.flow) then
  dfsw = one
else
  dfsw = f30*(one-t)*(t-one)*(t-one+eta)*(t-one+eta)/(eta**5)
endif
return
end function dfsw
!
subroutine ptcart(label,ncol,icol,x)
implicit none
!
! dump an array (ngrid,ncol) or just a column.
!
character (len=*), intent(in) :: label
integer, intent(in)           :: ncol, icol
real*8, dimension(ngrid,ncol), intent(in) :: x
!
integer :: ig, noff, nprt, ic, j
!
! print an header:
!
if (ncol.eq.1) then
  write (iout,'(3x,a,1x,"(column ",i4")")') label, icol
else
  write (iout,'(3x,a)') label
end if
if (ncol.eq.1) then
  do ig = 1, ngrid
    write(iout,1000) ig, x(ig,1)
  end do
!
else
  noff = mod (ncol,5)
  nprt = max(ncol - noff,0)
  do ic = 1, nprt, 5
    write(iout,1010) (j, j = ic, ic+4)
    do ig = 1, ngrid
      write(iout,1020) ig, x(ig,ic:ic+4)
    end do
  end do
  write (iout,1010) (j, j = nprt+1, nprt+noff)
  do ig = 1, ngrid
    write(iout,1020) ig, x(ig,nprt+1:nprt+noff)
  end do
end if
!
1000 format(1x,i5,f14.8)
1010 format(6x,5i14)
1020 format(1x,i5,5f14.8)
return
end subroutine ptcart
!
subroutine prtsph(label,ncol,icol,x)
implicit none
!
! dump an array (nbasis,ncol) or just a column.
!
character (len=*), intent(in) :: label
integer, intent(in)           :: ncol, icol
real*8, dimension(nbasis,ncol), intent(in) :: x
!
integer :: l, m, ind, noff, nprt, ic, j
!
! print an header:
!
if (ncol.eq.1) then
  write (iout,'(3x,a,1x,"(column ",i4")")') label, icol
else
  write (iout,'(3x,a)') label
end if
if (ncol.eq.1) then
  do l = 0, lmax
    ind = l*l + l + 1
    do m = -l, l
      write(iout,1000) l, m, x(ind+m,1)
    end do
  end do
!
else
  noff = mod (ncol,5)
  nprt = max(ncol - noff,0)
  do ic = 1, nprt, 5
    write(iout,1010) (j, j = ic, ic+4)
    do l = 0, lmax
      ind = l*l + l + 1
      do m = -l, l
        write(iout,1020) l, m, x(ind+m,ic:ic+4)
      end do
    end do
  end do
  write (iout,1010) (j, j = nprt+1, nprt+noff)
  do l = 0, lmax
    ind = l*l + l + 1
    do m = -l, l
      write(iout,1020) l, m, x(ind+m,nprt+1:nprt+noff)
    end do
  end do
end if
!
1000 format(1x,i3,i4,f14.8)
1010 format(8x,5i14)
1020 format(1x,i3,i4,5f14.8)
return
end subroutine prtsph
!
subroutine calcv(first,isph,g,pot,sigma,basloc,vplm,vcos,vsin)
implicit none
!
! computes a line of the (sparse) matrix/vector product for ddCOSMO.
!
logical, intent(in) :: first
integer, intent(in) :: isph
real*8, dimension(ngrid),       intent(in)    :: g
real*8, dimension(nbasis,nsph), intent(in)    :: sigma
real*8, dimension(ngrid),       intent(inout) :: pot
real*8, dimension(nbasis),      intent(inout) :: basloc, vplm
real*8, dimension(lmax+1),      intent(inout) :: vcos, vsin
!
integer :: ig, ij, jsph
real*8  :: vij(3), sij(3)
real*8  :: vvij, tij, xij, oij
!
pot = g
if (first) return
!
do ig = 1, ngrid
  if (ui(ig,isph).lt.one) then
    do ij = inl(isph), inl(isph+1) - 1
      jsph = nl(ij)
      vij  = csph(:,isph) + rsph(isph)*grid(:,ig) - csph(:,jsph)
      vvij = sqrt(dot_product(vij,vij))
      tij  = vvij/rsph(jsph) 
      if (tij.lt.one) then
        sij  = vij/vvij
        xij  = fsw(tij,eta*rsph(jsph))
        if (fi(ig,isph).gt.one) then
          oij = xij/fi(ig,isph)
        else
          oij = xij
        end if
        call ylmbas(sij,basloc,vplm,vcos,vsin)
        pot(ig) = pot(ig) + oij*intmlp(tij,sigma(:,jsph),basloc)
      end if
    end do
  end if
end do
!
return
end subroutine calcv
!
subroutine intrhs(isph,x,xlm)
implicit none
integer, intent(in) :: isph
real*8, dimension(ngrid),  intent(in)    :: x
real*8, dimension(nbasis), intent(inout) :: xlm
!
integer ig
xlm = zero
do ig = 1, ngrid
  xlm = xlm + basis(:,ig)*w(ig)*x(ig)
end do
!
if (iprint.ge.5) then
  call ptcart('pot',1,isph,x)
  call prtsph('vlm',1,isph,xlm)
end if
return
end subroutine intrhs
!
subroutine solve(isph,vlm,slm)
implicit none
integer, intent(in) :: isph
real*8, dimension(nbasis), intent(in)    :: vlm
real*8, dimension(nbasis), intent(inout) :: slm
!
slm = facl*vlm
!
if (iprint.ge.4) call prtsph('slm',1,isph,slm)
return
end subroutine solve
!
subroutine diis(n,nmat,x,e,b,xnew)
implicit none
integer,                             intent(in)    :: n
integer,                             intent(inout) :: nmat
real*8,  dimension(n,ndiis),         intent(inout) :: x, e
real*8,  dimension(ndiis+1,ndiis+1), intent(inout) :: b
real*8,  dimension(n),               intent(inout) :: xnew
!
integer :: nmat1, i
integer :: j, k
logical :: ok
!
real*8, allocatable :: bloc(:,:), cex(:)
!
if (nmat.ge.ndiis) then
  do j = 2, nmat - 10
    do k = 2, nmat - 10
      b(j,k) = b(j+10,k+10)
    end do
  end do
  do j = 1, nmat - 10
    x(:,j) = x(:,j+10)
    e(:,j) = e(:,j+10)
  end do
  nmat = nmat - 10
end if
nmat1 = nmat + 1
allocate (bloc(nmat1,nmat1),cex(nmat1))
memuse = memuse + (nmat1+1)*nmat1
memmax = max(memmax,memuse)
call makeb(n,nmat,e,b)
bloc   = b(1:nmat1,1:nmat1)
cex    = zero
cex(1) = one
call gjinv(nmat1,1,bloc,cex,ok)
if (.not. ok) then
  nmat = 1
  return
end if
xnew = zero
do i = 1, nmat
  xnew = xnew + cex(i+1)*x(:,i)
end do
nmat = nmat + 1
deallocate (bloc,cex)
memuse = memuse - (nmat1+1)*nmat1
return
end subroutine diis
!
subroutine makeb(n,nmat,e,b)
implicit none
integer, intent(in) :: n, nmat
real*8, dimension(n,ndiis),         intent(in) :: e
real*8, dimension(ndiis+1,ndiis+1), intent(inout) :: b
!
integer :: i
real*8  :: bij

if (nmat.eq.1) then
!
! first build:
!
  b(1,1) = zero
  b(1,2) = one
  b(2,1) = one
  b(2,2) = dot_product(e(:,1),e(:,1))
else
!
! first, update the lagrangian line:
!
  b(nmat+1,1) = one
  b(1,nmat+1) = one
!
! now, compute the new matrix elements:
!
  do i = 1, nmat - 1
    bij = dot_product(e(:,i),e(:,nmat))
    b(nmat+1,i+1) = bij
    b(i+1,nmat+1) = bij
  end do
  b(nmat+1,nmat+1) = dot_product(e(:,nmat),e(:,nmat))
end if
if (iprint.ge.5) then
  do i = 1, nmat + 1
    write(iout,'(21d12.4)') b(1:nmat+1,i)
  end do
end if
return
end subroutine makeb
!
subroutine ylmbas(x,basloc,vplm,vcos,vsin)
implicit none
real*8, dimension(3), intent(in) :: x
real*8, dimension(nbasis), intent(inout) :: basloc, vplm
real*8, dimension(lmax+1), intent(inout) :: vcos, vsin
!
integer :: l, m, ind
real*8  :: cthe, sthe, cphi, sphi, plm
!
! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
! coordinates of x.
!
cthe = x(3)
sthe = sqrt(one - cthe*cthe)
if (sthe.ne.zero) then
  cphi = x(1)/sthe
  sphi = x(2)/sthe
else
  cphi = zero
  sphi = zero
end if
!
! evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
! pointless if z = 1, as the only non vanishing terms will be the 
! ones with m=0.
!
if(sthe.ne.zero) then
  call trgev(cphi,sphi,vcos,vsin)
else
  vcos = one
  vsin = zero
end if
!
! evaluate the generalized legendre polynomials.
!
call polleg(cthe,sthe,vplm)
!
! now build the spherical harmonics. we will distinguish m=0,
! m>0 and m<0:
!
do l = 0, lmax
  ind = l**2 + l + 1
! m = 0
  basloc(ind) = facs(ind)*vplm(ind)
  do m = 1, l
    plm = vplm(ind+m)
!   m > 0
    basloc(ind+m) = facs(ind+m)*plm*vcos(m+1)
!   m < 0
    basloc(ind-m) = facs(ind-m)*plm*vsin(m+1)
  end do
end do
return
end subroutine ylmbas
!
subroutine dbasis(x,basloc,dbsloc,vplm,vcos,vsin)
implicit none
real*8, dimension(3),        intent(in)    :: x
real*8, dimension(nbasis),   intent(inout) :: basloc, vplm
real*8, dimension(3,nbasis), intent(inout) :: dbsloc
real*8, dimension(lmax+1),   intent(inout) :: vcos, vsin
!
integer :: l, m, ind
real*8  :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp
real*8  :: et(3), ep(3)
!
! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
! coordinates of x.
!
cthe = x(3)
sthe = sqrt(one - cthe*cthe)
if (sthe.ne.zero) then
  cphi = x(1)/sthe
  sphi = x(2)/sthe
else
  cphi = zero
  sphi = zero
end if
!
! evaluate the dirivatives of theta and phi:
!
et(1) = cthe*cphi
et(2) = cthe*sphi
et(3) = -sthe
if(sthe.ne.zero) then
  ep(1) = -sphi/sthe
  ep(2) = cphi/sthe
  ep(3) = zero
else
  ep(1) = zero
  ep(2) = zero
  ep(3) = zero
endif
!
! evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
! pointless if z = 1, as the only non vanishing terms will be the 
! ones with m=0.
!
if(sthe.ne.zero) then
  call trgev(cphi,sphi,vcos,vsin)
else
  vcos = one
  vsin = zero
end if
!
! evaluate the generalized legendre polynomials.
!
call polleg(cthe,sthe,vplm)
!
! now build the spherical harmonics. we will distinguish m=0,
! m>0 and m<0:
!
basloc = zero
dbsloc = zero
do l = 0, lmax
  ind = l*l + l + 1
! m = 0
  fln = facs(ind)
  basloc(ind) = fln*vplm(ind)
  if (l.gt.0) then
    dbsloc(:,ind) = fln*vplm(ind+1)*et(:)
  else
    dbsloc(:,ind) = zero
  end if
!dir$ simd
  do m = 1, l
    fln = facs(ind+m)
    plm = fln*vplm(ind+m)
    pp1 = zero
    if (m.lt.l) pp1 = -pt5*vplm(ind+m+1)
    pm1 = pt5*(dble(l+m)*dble(l-m+1)*vplm(ind+m-1))
    pp  = pp1 + pm1
!
!   m > 0
!
    basloc(ind+m)   = plm*vcos(m+1)
    dbsloc(:,ind+m) = -fln*pp*vcos(m+1)*et(:) - dble(m)*plm*vsin(m+1)*ep(:)
!
!   m < 0
!
    basloc(ind-m)   = plm*vsin(m+1)
    dbsloc(:,ind-m) = -fln*pp*vsin(m+1)*et(:) + dble(m)*plm*vcos(m+1)*ep(:)
  end do
end do
return
end subroutine dbasis
!
subroutine polleg(x,y,plm)
implicit none
real*8,                    intent(in)    :: x, y
real*8, dimension(nbasis), intent(inout) :: plm
!
! computes the l,m associated legendre polynomial for -1 <= x <= 1
! using the recurrence formula
!
!   (l-m)p(l,m) = x(2l-1)p(l-1,m) - (l+m-1)p(l-2,m)
!
integer :: m, ind, l, ind2
real*8  :: fact, pmm, somx2, pmm1, pmmo, pll, fm, fl
!
fact  = one
pmm   = one
somx2 = y
do m = 0, lmax 
  ind      = (m + 1)*(m + 1)
  plm(ind) = pmm
  if(m.eq.lmax) return
  fm = dble(m)
  pmm1 = x*(two*fm + one)*pmm
  ind2 = ind + 2*m + 2
  plm(ind2) = pmm1
  pmmo = pmm
  do l = m+2, lmax
    fl = dble(l)
    pll   = (x*(two*fl - one)*pmm1 - (fl + fm - one)*pmm)/(fl - fm)
    ind = l*l + l + 1 
    plm(ind+m) = pll
    pmm  = pmm1
    pmm1 = pll
  end do
  pmm  = -pmmo*fact*somx2
  fact = fact + two
end do
!
return
end subroutine polleg
!
subroutine trgev(x,y,cx,sx)
implicit none
real*8, intent(in) :: x, y
real*8, dimension(lmax+1), intent(inout) :: cx, sx
!
integer :: m
!
cx(1) = one
sx(1) = zero
cx(2) = x
sx(2) = y
do m = 3, lmax+1
  cx(m) = two*x*cx(m-1) - cx(m-2)
  sx(m) = two*x*sx(m-1) - sx(m-2)
end do
return
end subroutine trgev
!
real*8 function intmlp(t,sigma,basloc)
implicit none
real*8, intent(in) :: t
real*8, dimension(nbasis), intent(in) :: sigma, basloc
!
integer :: l, ind
real*8  :: tt, ss, fac
!
tt = one
ss = zero
do l = 0, lmax
  ind = l*l + l + 1
  fac = tt/facl(ind)
  ss = ss + fac*dot_product(basloc(ind-l:ind+l),sigma(ind-l:ind+l))
  tt = tt*t
end do
intmlp = ss
return
end function intmlp
!
subroutine itsolv(star,phi,psi,sigma,ene)
implicit none
logical,                        intent(in)    :: star
real*8, dimension(ncav),        intent(in)    :: phi
real*8, dimension(nbasis,nsph), intent(in)    :: psi
real*8,                         intent(inout) :: ene
real*8, dimension(nbasis,nsph), intent(inout) :: sigma
!
! local arrays:
!
real*8, allocatable :: g(:,:), pot(:), sigold(:,:), vlm(:), xi(:,:)
!
! scratch arrays:
!
real*8, allocatable :: basloc(:), vplm(:), vcos(:), vsin(:)
real*8, allocatable :: delta(:), norm(:)
!
! diis arrays:
!
real*8, allocatable :: xdiis(:,:,:), ediis(:,:,:), bmat(:)
!
! local variables:
!
integer :: it, isph, nmat, lenb, ig, c1, c2, cr
real*8  :: tol, drms, dmax, fep
logical :: dodiis, first
!
integer, parameter :: nitmax=300
real*8,  parameter :: ten=10.d0, tredis=1.0d-2
!
! initialize the timer:
!
call system_clock(count_rate=cr)
call system_clock(count=c1)
!
! allocate local variables and set convergence parameters
!
fep = (eps-one)/eps
tol = ten**(-iconv)
allocate (g(ngrid,nsph),pot(ngrid),vlm(nbasis),sigold(nbasis,nsph))
sigold = zero
allocate (delta(nbasis),norm(nsph))
allocate(vplm(nbasis),basloc(nbasis),vcos(lmax+1),vsin(lmax+1))
if (star) allocate(xi(ngrid,nsph))
memuse = memuse + ngrid*nsph + ngrid*nproc + 2*nbasis*nproc + nbasis*nsph + &
         2*nbasis*nproc + 2*(lmax+1)*nproc + nsph
if (star) memuse = memuse + nsph*ngrid
memmax = max(memmax,memuse)
!
! set up diis:
!
dodiis = .false.
nmat   = 1
lenb   = ndiis + 1
allocate (xdiis(nbasis,nsph,ndiis),ediis(nbasis,nsph,ndiis),bmat(lenb*lenb))
memuse = memuse + 2*nbasis*nsph*ndiis + lenb*lenb
memmax = max(memmax,memuse)
!
! solve the direct equations
!
if (.not. star) then
!
! build g:
!
  call wghpot(phi,g)
  do it = 1, nitmax
    first = it.eq.1
    vlm   = zero
!$omp parallel do default(shared) private(basloc,vcos,vsin,vplm) &
!$omp private(isph,pot,vlm,delta) schedule(dynamic,10) 
    do isph = 1, nsph
      call calcv(first,isph,g(:,isph),pot,sigold,basloc,vplm,vcos,vsin)
      call intrhs(isph,pot,vlm)
      call solve(isph,vlm,sigma(:,isph))
      delta = sigma(:,isph) - sigold(:,isph)
      call hsnorm(delta,norm(isph))
    end do
!$omp end parallel do
    call rmsvec(nsph,norm,drms,dmax)
    if (drms.le.tredis .and. ndiis.gt.0) dodiis = .true.
    if (dodiis) then
      xdiis(:,:,nmat) = sigma
      ediis(:,:,nmat) = sigma - sigold
      call diis(nbasis*nsph,nmat,xdiis,ediis,bmat,sigma)
    end if
    if (iprint.gt.1) then
      ene = pt5*fep*sprod(nbasis*nsph,sigma,psi)
      write(iout,1000) it, ene, drms, dmax
    end if
    if (drms.le.tol) goto 900
    sigold = sigma
  end do
  write(iout,1020) 
  stop
900 continue
  ene = pt5*fep*sprod(nbasis*nsph,sigma,psi)
!
else
!
! solve the adjoint equations:
!
  do it = 1, nitmax
    first = it.eq.1
    if (.not. first) then
      xi = zero
!$omp parallel do default(shared) private(isph,ig)
      do isph = 1, nsph
        do ig = 1, ngrid
          xi(ig,isph) = dot_product(sigold(:,isph),basis(:,ig))
        end do
      end do
!$omp end parallel do
    end if
!$omp parallel do default(shared) private(basloc,vcos,vsin,vplm) &
!$omp private(isph,vlm,delta) schedule(dynamic,10)
    do isph = 1, nsph
      call adjrhs(first,isph,psi(:,isph),xi,vlm,basloc,vplm,vcos,vsin)
      call solve(isph,vlm,sigma(:,isph))
      delta = sigma(:,isph) - sigold(:,isph)
      call hsnorm(delta,norm(isph))
    end do
!$omp end parallel do
    call rmsvec(nsph,norm,drms,dmax)
    if (drms.le.tredis .and. ndiis.gt.0) dodiis = .true.
    if (dodiis) then
      xdiis(:,:,nmat) = sigma
      ediis(:,:,nmat) = sigma - sigold
      call diis(nbasis*nsph,nmat,xdiis,ediis,bmat,sigma)
    end if
    if (iprint.gt.1) then
      write(iout,1000) it, zero, drms, dmax
    end if
    if (drms.le.tol) goto 910
    sigold = sigma
  end do
  write(iout,1020) 
  stop
910 continue
end if
!
! free the memory:
!
if (iprint.gt.1) write(iout,*)
if (star) deallocate(xi)
deallocate (g,pot,vlm,sigold)
deallocate (delta,norm)
deallocate (vplm,basloc,vcos,vsin)
deallocate (xdiis,ediis,bmat)
memuse = memuse - ngrid*nsph - ngrid*nproc - 2*nbasis*nproc - nbasis*nsph - &
         2*nbasis*nproc - 2*(lmax+1)*nproc - nsph
memuse = memuse - 2*nbasis*nsph*ndiis - lenb*lenb
!
call system_clock(count=c2)
if (iprint.gt.0) then
  if (star) then
    write(iout,1010) 'adjoint ', dble(c2-c1)/dble(cr)
  else
    write(iout,1010) '', dble(c2-c1)/dble(cr)
  end if
  write(iout,*)
end if
if (iprint.ge.3) then
  if (star) then
    call prtsph('solution to the ddcosmo adjoint equations:',nsph,0,sigma)
  else
    call prtsph('solution to the ddcosmo equations:',nsph,0,sigma)
  end if
end if
1000 format(' energy at iteration ',i4,': ',f14.7,' error (rms,max):',2f14.7)
1010 format(' the solution to the ddCOSMO ',a,'equations took ',f8.3,' seconds.')
1020 format(' ddCOSMO did not converge! Aborting...')
return
end subroutine itsolv
!
subroutine wghpot(phi,g)
implicit none
!
real*8, dimension(ncav),       intent(in)    :: phi
real*8, dimension(ngrid,nsph), intent(inout) :: g
!
integer isph, ig, ic
!
ic = 0
do isph = 1, nsph
  do ig = 1, ngrid
    if (ui(ig,isph).ne.zero) then
      ic = ic + 1
      g(ig,isph) = - ui(ig,isph)*phi(ic)
    end if
  end do
end do
!
return
end subroutine wghpot
!
subroutine hsnorm(u,unorm)
implicit none
real*8, dimension(nbasis), intent(in)    :: u
real*8,                    intent(inout) :: unorm
!
integer :: l, m, ind
real*8  :: fac
!
! compute the energy norm of a vector
!
unorm = zero
do l = 0, lmax
  ind = l*l + l + 1
  fac = one/(one + dble(l))
  do m = -l, l
    unorm = unorm + fac*u(ind+m)*u(ind+m)
  end do
end do
unorm = sqrt(unorm)
!
return
end subroutine hsnorm
!
subroutine adjrhs(first,isph,psi,xi,vlm,basloc,vplm,vcos,vsin)
implicit none
!
! compute a line of the (sparse) adjoint ddcosmo matrix/vector product.
!
logical,                       intent(in)    :: first
integer,                       intent(in)    :: isph
real*8, dimension(nbasis),     intent(in)    :: psi
real*8, dimension(ngrid,nsph), intent(in)    :: xi
real*8, dimension(nbasis),     intent(inout) :: vlm
real*8, dimension(nbasis),     intent(inout) :: basloc, vplm
real*8, dimension(lmax+1),     intent(inout) :: vcos, vsin
!
integer :: ij, jsph, ig, l, ind, m
real*8  :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t
vlm = psi
if (first) return
!
do ij = inl(isph), inl(isph+1) - 1
  jsph = nl(ij)
  do ig = 1, ngrid
    vji  = csph(:,jsph) + rsph(jsph)*grid(:,ig) - csph(:,isph)
    vvji = sqrt(dot_product(vji,vji))
    tji  = vvji/rsph(isph)
    if (tji.lt.one) then
      sji = vji/vvji
      xji = fsw(tji,eta*rsph(isph))
      if (fi(ig,jsph).gt.one) then
        oji = xji/fi(ig,jsph)
      else
        oji = xji
      end if
      call ylmbas(sji,basloc,vplm,vcos,vsin)
      t   = one
      fac = w(ig)*xi(ig,jsph)*oji
      do l = 0, lmax
        ind  = l*l + l + 1
        ffac = fac*t/facl(ind)
        do m = -l, l
          vlm(ind+m) = vlm(ind+m) + ffac*basloc(ind+m)
        end do
        t = t*tji
      end do
    end if
  end do
end do
!
return
end subroutine adjrhs
!
subroutine rmsvec(n,v,vrms,vmax)
implicit none
integer,               intent(in)    :: n
real*8,  dimension(n), intent(in)    :: v
real*8,                intent(inout) :: vrms, vmax
!
integer i
vrms = zero
vmax = zero
do i = 1, n
  vmax = max(vmax,abs(v(i)))
  vrms = vrms + v(i)*v(i)
end do
vrms = sqrt(vrms/dble(n))
return
end subroutine rmsvec
!
subroutine gjinv(n,nrhs,a,b,ok)
implicit none
!
integer,                    intent(in)    :: n, nrhs
logical,                    intent(inout) :: ok
real*8,  dimension(n,n),    intent(inout) :: a
real*8,  dimension(n,nrhs), intent(inout) :: b
!
integer :: i, j, k, irow, icol
real*8  :: big, dum, pinv
!
integer, allocatable :: indxc(:), indxr(:), piv(:)
real*8,  allocatable :: scr(:)
!
allocate (indxc(n), indxr(n), piv(n))
allocate (scr(n))
!
memuse = memuse + 4*n
memmax = max(memmax,memuse)
!
ok  = .false.
piv = 0
!
irow = 0
icol = 0
do i = 1, n
  big = zero
  do j = 1, n
    if (piv(j).ne.1) then
      do k = 1, n
        if (piv(k).eq.0) then
          if (abs(a(j,k)).gt.big) then
            big  = abs(a(j,k))
            irow = j
            icol = k
          end if
        end if
      end do
    end if
  end do
!
  piv(icol) = piv(icol) + 1
  if (piv(icol) .gt. 1) then
    write(iout,1000)
    return
  end if
  if (irow.ne.icol) then
    scr         = a(irow,:)
    a(irow,:)   = a(icol,:)
    a(icol,:)   = scr  
    scr(1:nrhs) = b(irow,:)
    b(irow,:)   = b(icol,:)
    b(icol,:)   = scr(1:nrhs)       
  end if
!
  indxr(i) = irow
  indxc(i) = icol
!
  if (a(icol,icol) .eq. zero) then
    write(iout,1000)
    return
  end if
!
  pinv = one/a(icol,icol)
  a(icol,icol) = one
  a(icol,:) = a(icol,:)*pinv
  b(icol,:) = b(icol,:)*pinv
!
  do j = 1, n
    if (j.ne.icol) then
      dum       = a(j,icol)
      a(j,icol) = zero
      a(j,:)    = a(j,:) - a(icol,:)*dum
      b(j,:)    = b(j,:) - b(icol,:)*dum
    end if
  end do
end do
!
do j = n, 1, -1
  if (indxr(j) .ne. indxc(j)) then
    scr           = a(:,indxr(j))
    a(:,indxr(j)) = a(:,indxc(j))
    a(:,indxc(j)) = scr
  end if
end do
!
ok = .true.
deallocate (indxr,indxc,piv,scr)
memuse = memuse - 4*n
return
!
1000 format (' warning: singular matrix in gjinv!')
end subroutine gjinv
!
subroutine header
implicit none
!
1000 format( /,&
             '      888      888  .d8888b.   .d88888b.   .d8888b.  888b     d888  .d88888b. ',/,  &
             '      888      888 d88P  Y88b d88P" "Y88b d88P  Y88b 8888b   d8888 d88P" "Y88b',/,  &
             '      888      888 888    888 888     888 Y88b.      88888b.d88888 888     888',/,  &
             '  .d88888  .d88888 888        888     888  "Y888b.   888Y88888P888 888     888',/,  &
             ' d88" 888 d88" 888 888        888     888     "Y88b. 888 Y888P 888 888     888',/,  &
             ' 888  888 888  888 888    888 888     888       "888 888  Y8P  888 888     888',/,  &
             ' Y88b 888 Y88b 888 Y88b  d88P Y88b. .d88P Y88b  d88P 888   "   888 Y88b. .d88P',/,  &
             '  "Y88888  "Y88888  "Y8888P"   "Y88888P"   "Y8888P"  888       888  "Y88888P" ',/,  &
             '                                                                              ',/,  &
             ' An implementation of COSMO using a domain decomposition linear scaling strategy.',/)
1010 format( ' parameters:',/, &
             '   number of grid points:                  ',i8,/,   &
             '   number of spheres:                      ',i8,/,   &
             '   lmax for the spherical harmonics basis: ',i8,/,   &
             '   convergence threshold:                  ',d8.1,/, &
             '   regularization parameter:               ',f8.2,/)
if (iprint.gt.0) write(iout,1000)
if (iprint.gt.0) write(iout,1010) ngrid, nsph, lmax, 10.0d0**(-iconv), eta
return
end subroutine header
!
subroutine fdoka(isph,sigma,xi,basloc,dbsloc,vplm,vcos,vsin,fx)
implicit none
integer,                         intent(in)    :: isph
real*8,  dimension(nbasis,nsph), intent(in)    :: sigma
real*8,  dimension(ngrid),       intent(in)    :: xi
real*8,  dimension(nbasis),      intent(inout) :: basloc, vplm
real*8,  dimension(3,nbasis),    intent(inout) :: dbsloc
real*8,  dimension(lmax+1),      intent(inout) :: vcos, vsin
real*8,  dimension(3),           intent(inout) :: fx
!
integer :: ig, ij, jsph, l, ind, m
real*8  :: vvij, tij, xij, oij, t, fac, fl, f1, f2, f3, beta
real*8  :: vij(3), sij(3), alp(3), va(3)
!
do ig = 1, ngrid
  va = zero
  do ij = inl(isph), inl(isph+1) - 1
    jsph = nl(ij)
    vij  = csph(:,isph) + rsph(isph)*grid(:,ig) - csph(:,jsph)
    vvij = sqrt(dot_product(vij,vij))
    tij  = vvij/rsph(jsph)
    if (tij.ge.one) cycle
    sij  = vij/vvij
    call dbasis(sij,basloc,dbsloc,vplm,vcos,vsin)
    alp  = zero
    t    = one
    do l = 1, lmax
      ind = l*l + l + 1
      fl  = dble(l)
      fac = t/facl(ind)
      do m = -l, l
        f2 = fac*sigma(ind+m,jsph)
        f1 = f2*fl*basloc(ind+m)
        alp(:) = alp(:) + f1*sij(:) + f2*dbsloc(:,ind+m)
      end do
      t = t*tij
    end do
    beta = intmlp(tij,sigma(:,jsph),basloc)
    xij = fsw(tij,eta*rsph(jsph))
    if (fi(ig,isph).gt.one) then
      oij = xij/fi(ig,isph)
      f2  = -oij/fi(ig,isph)
    else
      oij = xij
      f2  = zero
    end if
    f1 = oij/rsph(jsph)
    va(:) = va(:) + f1*alp(:) + beta*f2*zi(:,ig,isph)
    if (tij .gt. (one-eta*rsph(jsph))) then
      f3 = beta*dfsw(tij,eta*rsph(jsph))/rsph(jsph)
      if (fi(ig,isph).gt.one) f3 = f3/fi(ig,isph)
      va(:) = va(:) + f3*sij(:)
    end if
  end do
  fx = fx - w(ig)*xi(ig)*va(:)
end do
return
end subroutine fdoka
!
subroutine fdokb(isph,sigma,xi,basloc,dbsloc,vplm,vcos,vsin,fx)
implicit none
integer,                         intent(in)    :: isph
real*8,  dimension(nbasis,nsph), intent(in)    :: sigma
real*8,  dimension(ngrid,nsph),  intent(in)    :: xi
real*8,  dimension(nbasis),      intent(inout) :: basloc, vplm
real*8,  dimension(3,nbasis),    intent(inout) :: dbsloc
real*8,  dimension(lmax+1),      intent(inout) :: vcos, vsin
real*8,  dimension(3),           intent(inout) :: fx
!
integer :: ig, ji, jsph, l, ind, m, jk, ksph
logical :: proc
real*8  :: vvji, tji, xji, oji, t, fac, fl, f1, f2, beta, di
real*8  :: b, g1, g2, vvjk, tjk, f, xjk
real*8  :: vji(3), sji(3), alp(3), vb(3), vjk(3), sjk(3), vc(3)
!
do ig = 1, ngrid
  vb = zero
  vc = zero
  do ji = inl(isph), inl(isph+1) - 1
    jsph = nl(ji)
    vji  = csph(:,jsph) + rsph(jsph)*grid(:,ig) - csph(:,isph)
    vvji = sqrt(dot_product(vji,vji))
    tji  = vvji/rsph(isph)
    if (tji.gt.one) cycle
    sji  = vji/vvji
    call dbasis(sji,basloc,dbsloc,vplm,vcos,vsin)
!
    alp = zero
    t   = one
    do l = 1, lmax
      ind = l*l + l + 1
      fl  = dble(l)
      fac = t/facl(ind)
      do m = -l, l
        f2 = fac*sigma(ind+m,isph)
        f1 = f2*fl*basloc(ind+m)
        alp = alp + f1*sji + f2*dbsloc(:,ind+m)
      end do
      t = t*tji
    end do
    xji = fsw(tji,eta*rsph(isph))
    if (fi(ig,jsph).gt.one) then
      oji = xji/fi(ig,jsph)
    else
      oji = xji
    end if
    f1 = oji/rsph(isph)
    vb = vb + f1*alp*xi(ig,jsph)
    if (tji.gt.one-eta*rsph(isph)) then
      beta = intmlp(tji,sigma(:,isph),basloc)
      if (fi(ig,jsph) .gt. one) then
        di  = one/fi(ig,jsph)
        fac = di*xji
        proc = .false.
        b    = zero
        do jk = inl(jsph), inl(jsph+1) - 1
          ksph = nl(jk)
          vjk  = csph(:,jsph) + rsph(jsph)*grid(:,ig) - csph(:,ksph)
          vvjk = sqrt(dot_product(vjk,vjk))
          tjk  = vvjk/rsph(ksph)
          if (ksph.ne.isph) then
            if (tjk.lt.one) then
              proc = .true.
              sjk  = vjk/vvjk
              call ylmbas(sjk,basloc,vplm,vcos,vsin)
              g1  = intmlp(tjk,sigma(:,ksph),basloc)
              xjk = fsw(tjk,eta*rsph(ksph))
              b   = b + g1*xjk
            end if
          end if
        end do
        if (proc) then
          g1 = di*di*dfsw(tji,eta*rsph(isph))/rsph(isph)
          g2 = g1*xi(ig,jsph)*b
          vc = vc + g2*sji
        end if
      else
        di  = one
        fac = zero
      end if
      f2 = (one-fac)*di*dfsw(tji,eta*rsph(isph))/rsph(isph)
      vb = vb + f2*xi(ig,jsph)*beta*sji
    end if 
  end do
  fx = fx + w(ig)*(vb - vc)
! fx = fx - w(ig)*vc
end do
return
end subroutine fdokb
!
subroutine fdoga(isph,xi,phi,fx)
implicit none
integer,                        intent(in)    :: isph
real*8,  dimension(ngrid,nsph), intent(in)    :: xi, phi
real*8,  dimension(3),          intent(inout) :: fx
!
integer :: ig, ji, jsph
real*8  :: vvji, tji, fac, swthr
real*8  :: alp(3), vji(3), sji(3)
!
do ig = 1, ngrid
  alp = zero
  if (ui(ig,isph) .gt. zero .and. ui(ig,isph).lt.one) then
    alp = alp + phi(ig,isph)*xi(ig,isph)*zi(:,ig,isph)
  end if
  do ji = inl(isph), inl(isph+1) - 1
    jsph  = nl(ji)
    vji   = csph(:,jsph) + rsph(jsph)*grid(:,ig) - csph(:,isph)
    vvji  = sqrt(dot_product(vji,vji))
    tji   = vvji/rsph(isph)
    swthr = one - eta*rsph(isph)
    if (tji.lt.one .and. tji.gt.swthr .and. ui(ig,jsph).gt.zero) then
      sji = vji/vvji
      fac = - dfsw(tji,eta*rsph(isph))/rsph(isph)
      alp = alp + fac*phi(ig,jsph)*xi(ig,jsph)*sji
    end if
  end do
  fx = fx - w(ig)*alp
end do
return 
end subroutine fdoga
!
end module ddcosmo
