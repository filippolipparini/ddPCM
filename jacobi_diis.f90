!---------------------------------------------------------------------------------------------
! Purpose : Jacobi/DIIS solver
!
! variables:
!
!   n        : integer, input, size of the matrix
!
!   lprint   : integer, input, printing flag.
!
!   diis_max : integer, input, number of points to be used for diis extrapolation
!
!              if diis_max = 0, this is just a Jacobi solver.
!
!   norm     : integer, input, norm to be used to evaluate convergence
!              1: max |x_new - x|
!              2: rms (x_new - x)
!              3: rms (x_new - x) and max |x_new - x|
!              4: norm computed by the user-provided function u_norm(n,x)
!
!   tol      : real, input, convergence criterion. if norm = 3, convergence is 
!              achieved when rms (x_new - x) < tol and max |x_new - x| < 10*tol.
!
!   rhs      : real, dimension(n), input, right-hand side of the linear system
!
!   x        : real, dimension(n). In input, a guess of the solution (can be zero).
!              In output, the solution
!
!   n_iter   : integer, in input, the maximum number of iterations. In output,
!              the number of iterations needed to converge.
!
!   ok       : logical, output, T if the solver converged, false otherwise.
!
!   matvec   : external, subroutine to compute the required matrix-vector multiplication
!              format: subroutine matvec(n,x,y)
!
!   dm1vec   : external, subroutine to apply the inverse diagonal matrix to a vector.
!              format: subroutine dm1vec(n,x,y)
! 
!   u_norm   : external, optional function to compute the norm of a vector.
!              format: real*8 function u_norm(n,x)
!
!---------------------------------------------------------------------------------------------
subroutine jacobi_diis( n, lprint, diis_max, norm, tol, rhs, x, n_iter, ok, matvec, &
                        dm1vec, u_norm)
!                
      implicit none
      integer,               intent(in)    :: n, diis_max, norm, lprint
      real*8,                intent(in)    :: tol
      real*8,  dimension(n), intent(in)    :: rhs
      real*8,  dimension(n), intent(inout) :: x
      integer,               intent(inout) :: n_iter
      logical,               intent(inout) :: ok
      external                             :: matvec, dm1vec
      real*8,  optional                    :: u_norm
      external                             :: u_norm
!
      integer :: it, nmat, istatus, lenb
      real*8  :: rms_norm, max_norm, tol_max
      logical :: dodiis
!
      real*8, allocatable :: x_new(:), y(:), x_diis(:,:), e_diis(:,:), bmat(:,:)
!
!---------------------------------------------------------------------------------------------
!
!     check inputs
      if ( (norm.eq.4) .and. (.not.present(u_norm)) ) then
        write(6,*) ' must provide a function norm(n,x) to evaluate the norm of the increment'
        stop
      endif
!
!     DIIS extrapolation flag
      dodiis =  (diis_max.ne.0)
!
!     set tolerance
      tol_max = 10.0d0 * tol
!
!     extrapolation required
      if (dodiis) then
!
!       allocate workspaces
        lenb = diis_max + 1
        allocate( x_diis(n,diis_max), e_diis(n,diis_max), bmat(lenb,lenb) , stat=istatus )
        if (istatus .ne. 0) then
          write(*,*) ' jacobi_diis: [1] failed allocation (diis)'
          stop
        endif
!        
!       an enigmatic constant
        nmat = 1
!        
      endif
!
!     allocate workspaces
      allocate( x_new(n), y(n) , stat=istatus )
      if (istatus .ne. 0) then
        write(*,*) ' jacobi_diis: [2] failed allocation (scratch)' 
        stop
      endif
!
!     Jacobi iterations
!     =================
      do it = 1, n_iter
!
!       y = rhs - O x
        call matvec( n, x, y )
        y = rhs - y
!
!       x_new = D^-1 y
        call dm1vec(n,y,x_new)
!
!       DIIS extrapolation
!       ==================
        if (dodiis) then
!
          x_diis(:,nmat) = x_new
          e_diis(:,nmat) = x_new - x
!
          call diis(n,nmat,diis_max,x_diis,e_diis,bmat,x_new)
!
        endif
!
!       increment
        x = x_new - x
!
!       rms/max norm of increment
        if ( norm.le.3 ) then
!
!         compute norm
          call rmsvec( n, x, rms_norm, max_norm )
!
!         check norm
          if ( norm.eq.1 ) then
!                  
            ok = (rms_norm.lt.tol)
            
          elseif ( norm.eq.2 ) then
!                  
            ok = (max_norm.lt.tol)
!            
          else 

            ok = (rms_norm.lt.tol) .and. (max_norm.lt.tol_max)
!            
          endif
!
!       user-provided norm of increment
        elseif ( norm.eq.4 ) then
!
!         just a placeholder for printing
          max_norm = -1.d0
!
!         compute norm
          rms_norm = u_norm( n, x )
!
!         check norm
          ok = (rms_norm.lt.tol)
!          
        endif
!
!       printing
        if ( lprint.gt.0 ) then
           if (norm.eq.1) then
             write(*,110) it, 'max', max_norm
           else if (norm.eq.2) then
             write(*,110) it, 'rms', rms_norm
           else if (norm.eq.3) then
             write(*,100) it, rms_norm, max_norm
           else if (norm.eq.4) then
             write(*,120) it, rms_norm
           end if
         end if
  100   format(t3,'iter=',i4,' residual norm (rms,max): ', 2d14.4 )
  110   format(t3,'iter=',i4,' residual norm (',a,'): ', d14.4 )
  120   format(t3,'iter=',i4,' residual norm: ', d14.4 )
!
!       update
        x = x_new
!
!       EXIT Jacobi loop here
!       =====================
        if (ok) exit
!
      enddo
!
!     record number of Jacobi iterations
      n_iter = it
!
      return
!
!
endsubroutine jacobi_diis
!---------------------------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------------------------
!
subroutine diis(n,nmat,ndiis,x,e,b,xnew)
implicit none
integer,                             intent(in)    :: n, ndiis
integer,                             intent(inout) :: nmat
real*8,  dimension(n,ndiis),         intent(inout) :: x, e
real*8,  dimension(ndiis+1,ndiis+1), intent(inout) :: b
real*8,  dimension(n),               intent(inout) :: xnew
!
integer :: nmat1, i, istatus
integer :: j, k
logical :: ok
!
real*8, allocatable :: bloc(:,:), cex(:)
!
real*8, parameter :: zero = 0.0d0, one = 1.0d0
!
!------------------------------------------------------------------------------
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
allocate (bloc(nmat1,nmat1),cex(nmat1) , stat=istatus)
if ( istatus.ne.0 ) then 
  write(*,*) 'diis: allocation failed!'
  stop
endif

call makeb(n,nmat,ndiis,e,b)
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
deallocate (bloc,cex , stat=istatus)
if ( istatus.ne.0 ) then 
  write(*,*) 'diis: deallocation failed!'
  stop
endif
!
return
end subroutine diis
  !
subroutine makeb(n,nmat,ndiis,e,b)
implicit none
integer, intent(in) :: n, nmat, ndiis
real*8, dimension(n,ndiis),         intent(in) :: e
real*8, dimension(ndiis+1,ndiis+1), intent(inout) :: b
!
integer :: i
real*8  :: bij
real*8, parameter :: zero = 0.0d0, one = 1.0d0
  
! 1st built
if (nmat.eq.1) then
!
!       [ 0 |  1  ]
!   b = [ --+---- ]
!       [ 1 | e*e ]
!
  b(1,1) = zero
  b(1,2) = one
  b(2,1) = one
  b(2,2) = dot_product(e(:,1),e(:,1))
!
! subsequent builts
else
!
!   first, update the lagrangian line:
  b(nmat+1,1) = one
  b(1,nmat+1) = one
!
!   now, compute the new matrix elements:
  do i = 1, nmat - 1
    bij = dot_product(e(:,i),e(:,nmat))
    b(nmat+1,i+1) = bij
    b(i+1,nmat+1) = bij
  end do
  b(nmat+1,nmat+1) = dot_product(e(:,nmat),e(:,nmat))
end if
!
return
end subroutine makeb
!
subroutine gjinv(n,nrhs,a,b,ok)
implicit none
!
integer,                    intent(in)    :: n, nrhs
logical,                    intent(inout) :: ok
real*8,  dimension(n,n),    intent(inout) :: a
real*8,  dimension(n,nrhs), intent(inout) :: b
!
integer :: i, j, k, irow, icol, istatus
real*8  :: big, dum, pinv
real*8, parameter :: zero = 0.0d0, one = 1.0d0
!
integer, allocatable :: indxc(:), indxr(:), piv(:)
real*8,  allocatable :: scr(:)
!
allocate (indxc(n), indxr(n), piv(n) , stat=istatus)
if ( istatus.ne.0 ) then
  write(*,*)'gjinv: allocation failed! [1]'
  stop
endif
allocate (scr(n) , stat=istatus)
if ( istatus.ne.0 ) then
  write(*,*)'gjinv: allocation failed! [2]'
  stop
endif
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
    write(*,1000)
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
    write(*,1000)
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
deallocate (indxr,indxc,piv,scr , stat=istatus)
if ( istatus.ne.0 ) then
  write(*,*)'gjinv: deallocation failed! [1]'
  stop
endif

return
1000 format (' warning: singular matrix in gjinv!')
end subroutine gjinv
!
!------------------------------------------------------------------------------
! Purpose : compute root-mean-square and max norm
!------------------------------------------------------------------------------
subroutine rmsvec( n, v, vrms, vmax )
!
      implicit none
      integer,               intent(in)    :: n
      real*8,  dimension(n), intent(in)    :: v
      real*8,                intent(inout) :: vrms, vmax
!
      integer :: i
      real*8, parameter :: zero=0.0d0
!      
!------------------------------------------------------------------------------
!      
!     initialize
      vrms = zero
      vmax = zero
!
!     loop over entries
      do i = 1,n
!
!       max norm
        vmax = max(vmax,abs(v(i)))
!
!       rms norm
        vrms = vrms + v(i)*v(i)
!        
      enddo
!
!     the much neglected square root
      vrms = sqrt(vrms/dble(n))
!      
      return
!      
!      
endsubroutine rmsvec
!------------------------------------------------------------------------------
