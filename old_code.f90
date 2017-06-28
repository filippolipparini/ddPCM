!-------------------------------------------------------------------------------
! Purpose : solve COSMO equation 
!
!             L sigma = g
!
!           or adjoint equation
!
!             L^* sigma = Psi
!
!           Compute solvation energy
!
!             E_s = 1/2 f(eps) < Psi , sigma >
!
!-------------------------------------------------------------------------------
subroutine old_itsolv( star, phi, psi, sigma, ene )
      use ddcosmo
!
      implicit none
      logical,                        intent(in)    :: star
      real*8, dimension(ncav),        intent(in)    :: phi
      real*8, dimension(nbasis,nsph), intent(in)    :: psi
      real*8,                         intent(inout) :: ene
      real*8, dimension(nbasis,nsph), intent(inout) :: sigma
!
!     local arrays:
      real*8, allocatable :: g(:,:), pot(:), sigold(:,:), vlm(:), xi(:,:)
!
!     scratch arrays:
      real*8, allocatable :: basloc(:), vplm(:), vcos(:), vsin(:)
      real*8, allocatable :: delta(:), norm(:)
!
!     diis arrays:
      real*8, allocatable :: xdiis(:,:,:), ediis(:,:,:), bmat(:)
!
!     local variables:
      integer :: it, isph, nmat, lenb, ig, c1, c2, cr, j
      real*8  :: tol, drms, dmax, fep, s1, s2, s3
      logical :: dodiis, first
!
      integer, parameter :: nitmax=300
      real*8,  parameter :: ten=10.d0, tredis=1.0d-2
      integer :: istatus
!
      real*8 :: Lsigma(nbasis), gvoid(nbasis,nsph),drmsold
!      
!-------------------------------------------------------------------------------
!
!     initialize the timer:
      call system_clock(count_rate=cr)
      call system_clock(count=c1)
!
!     allocate local variables and set convergence parameters
!
!     compute f(eps) 
      fep = (eps-one)/eps
!
!     threshold for solver
      tol = ten**(-iconv)
!
      allocate( g(ngrid,nsph), pot(ngrid), vlm(nbasis), sigold(nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [1] failed allocation !'
        stop
      endif
!  
!     initialize
      sigold(:,:)=zero
!  
      allocate( delta(nbasis), norm(nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [2] failed allocation !'
        stop
      endif
!
!     initialize
      norm(:)=zero
!
      allocate( vplm(nbasis), basloc(nbasis), vcos(lmax+1), vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [3] failed allocation !'
        stop
      endif
!
!     adjoint equation
      if (star) then
        allocate( xi(ngrid,nsph) , stat=istatus )
        if ( istatus.ne.0 ) then
          write(*,*)'itsolv : [4] failed allocation !'
          stop
        endif
      endif
!
!     update memory usage
      memuse = memuse + ngrid*nsph + ngrid*nproc + 2*nbasis*nproc + nbasis*nsph + &
               2*nbasis*nproc + 2*(lmax+1)*nproc + nsph
      if (star) memuse = memuse + nsph*ngrid
      memmax = max(memmax,memuse)
!
!     set up diis:
      dodiis = .false.
      nmat   = 1
      lenb   = ndiis + 1
      allocate( xdiis(nbasis,nsph,ndiis), ediis(nbasis,nsph,ndiis), bmat(lenb*lenb) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [5] failed allocation !'
        stop
      endif
!
!     update memory usage
      memuse = memuse + 2*nbasis*nsph*ndiis + lenb*lenb
      memmax = max(memmax,memuse)
!
!
!     solve the direct equation S \sigma = \Phi
!     -----------------------------------------
!
!     numerical discretization :
!
!     [ L_11 ... L_1M ] [ X_1 ]   [ g_1 ]
!     [  .        .   ] [  .  ]   [  .  ]
!     [  .        .   ] [  .  ] = [  .  ]
!     [  .        .   ] [  .  ]   [  .  ]
!     [ L_M1 ... L_MM ] [ X_M ]   [ g_M ]
!
!     or equivalently :
!
!     L_ii X_i = g_i - sum_j L_ij X_j
!
!     using Jacobi method :
!
!             n                       n-1
!     L_ii X_i  = g_i - sum_j L_ij X_j
!
      if (.not. star) then
!
!       print header
        if (iprint.gt.1) then
          write(iout,*)'Solution of DIRECT system:'      
        endif
!
!       build g
        g(:,:) = zero
        call wghpot( phi, g )
!
!       iterate of Jacobi method
        do it = 1, nitmax
!
!         1st iteration flag
          first = ( it.eq.1 )
!
!         initialize vlm = g_i - sum_j L_ij X_j
          vlm = zero
!    
      !$omp parallel do default(shared) private(basloc,vcos,vsin,vplm) &
      !$omp private(isph,pot,vlm,delta) schedule(dynamic,10) 
!
!         loop over spheres
          do isph = 1, nsph
!
!           compute (sparse) matrix/vector product :
!
!           pot(:) = S \sigma_i^n  = sum_j \omega_ij * S_ij^n-1 \sigma_j + \Phi(:)
!
            call calcv( first, isph, g(:,isph), pot, sigold, basloc, vplm, vcos, vsin )
!      
!           integrate pot(:) against spherical harmonics to obtain numerical discretization :
!
!           vlm(:) = g_i - sum_j L_ij X_j
!
            call intrhs( isph, pot, vlm )
!      
!           solve : L_ii X_i = vlm  ( L_ii is a diagonal matrix )
            call solve( isph, vlm, sigma(:,isph) )
!
!
!!!!===================================================================
!!!!           OLD VERSION
!!!!           compute error
!!!            delta(:) = sigma(:,isph) - sigold(:,isph)
!!!!
!!!!           energy norm of error      
!!!            call hsnorm( delta(:), norm(isph) )
!!!!===================================================================
!
!
          end do
!    
      !$omp end parallel do
!
!
!!!!===================================================================
!!!!         OLD VERSION
!!!!         compute root-mean-square and max norm
!!!          call rmsvec( nsph, norm, drms, dmax )
!!!!    
!!!          if ( drms.le.tredis .and. ndiis.gt.0 )  dodiis = .true.
!!!!===================================================================
!
!
!===================================================================
!         NEW VERSION
          if ( ndiis.gt.0 )  dodiis = .true.
!===================================================================
!          
!
          if ( dodiis ) then
!                  
            xdiis(:,:,nmat) = sigma
            ediis(:,:,nmat) = sigma - sigold
!            
            call diis(nbasis*nsph,nmat,xdiis,ediis,bmat,sigma)
!            
          endif
!
!
!===================================================================
!         NEW VERSION
          do isph = 1,nsph
!
!           compute error
            delta(:) = sigma(:,isph) - sigold(:,isph)
!
!           energy norm of error      
            call hsnorm( delta(:), norm(isph) )
!
          end do

          drmsold = drms
          call rmsvec( nsph, norm, drms, dmax )
!!!          if ( drms > drmsold)  write(iout,1030)
 1030     format(' residual does not decrease monotonically!')         
!===================================================================
!
!
!         printing
          if (iprint.gt.1) then
            ene = pt5*fep*sprod(nbasis*nsph,sigma,psi)
            write(iout,1000) it, ene, drms, dmax
          end if
!
!         tolerance has been reached, exit loop
          if (drms.le.tol) goto 900
!
!         tolerance has NOT been reached, update solution
          sigold = sigma
!
!       end of Jacobi iterations
        end do
!
!       UNSUCCESSFUL solve, end here
        write(iout,1020) 
        stop
!
!       SUCCESSFUL solve, just continue
  900   continue
!
!       compute energy :
!
!       E_s = 1/2 f(\epsilon_s) * < \Psi , \sigma >
!
        ene = pt5*fep*sprod(nbasis*nsph,sigma,psi)
!
! 
!       check solution
!       ==============
        if (iprint.gt.1) then
!
!         compute action
          s1 = zero ; s2 = zero ; s3 = zero ; pot = zero ; Lsigma = zero ; gvoid = zero
!          
          do isph = 1,nsph
!          
!           action of g_i - L_ij , j \ne i
            call calcv( first, isph, g(:,isph), pot, sigma, basloc, vplm, vcos, vsin )
            call intrhs( isph, pot, Lsigma )
            Lsigma = Lsigma
!            
!           action of L_ii
            Lsigma = Lsigma - sigma(:,isph) / facl
!          
            do j=1,nbasis
!            
              s1 = s1 + ( Lsigma(j) )**2
              s3 = s3 + ( sigma(j,isph) )**2
!              
            enddo
!
            do j=1,ngrid
              s2 = s2 + ( g(j,isph) )**2
            enddo
!
          enddo
!
          write(*,*) 'ddCOSMO : '
          write(*,*) ' '
          write(*,1002) sqrt(s3) , sqrt(s1) / sqrt(s2)
 1002     format(' || sigma || , || Phi - L sigma || / || Phi || =  ',2(e12.5,2x) )    
! 
          if ( abs(sqrt(s1)/sqrt(s2)) .gt. 1.E-10 ) then
            write(*,*) 'solution failed!'
            write(*,*)'PAUSE - type any key to continue'
            read(*,*)
          endif        
! 
        endif
!
!
!
!     solve the adjoint equations S' \sigma = \Psi
!     --------------------------------------------
!
      else
!
!       print header
        if (iprint.gt.1) then
          write(iout,*)'Solution of ADJOINT system:'     
        endif
!
!       Jacobi method iterations
        do it = 1, nitmax
!
!         1st iteration flag
          first = it.eq.1
!                             
!                                               n-1
!         update xi( int. point , j ) = \sigma_j    ( int. point )
          if (.not. first) then
!
            xi = zero
!          
      !$omp parallel do default(shared) private(isph,ig)
!          
            do isph = 1, nsph
              do ig = 1, ngrid
                xi(ig,isph) = dot_product(sigold(:,isph),basis(:,ig))
              end do
            end do
!          
      !$omp end parallel do
!          
          end if
!          
      !$omp parallel do default(shared) private(basloc,vcos,vsin,vplm) &
      !$omp private(isph,vlm,delta) schedule(dynamic,10)
!          
!         loop over spheres
          do isph = 1, nsph
!
!           compute rhs of adjoint problem
            call adjrhs(first,isph,psi(:,isph),xi,vlm,basloc,vplm,vcos,vsin)
!
!           solve L_ii X_i = vlm  ( matrix L_ii is diagonal )
            call solve(isph,vlm,sigma(:,isph))
!            
!
!!!!===================================================================
!!!!           OLD VERSION
!!!            delta = sigma(:,isph) - sigold(:,isph)
!!!            call hsnorm(delta,norm(isph))
!!!!===================================================================
!            
!            
          end do
!            
      !$omp end parallel do
!            
!!!!===================================================================
!!!!         OLD VERSION
!!!          call rmsvec(nsph,norm,drms,dmax)
!!!          if (drms.le.tredis .and. ndiis.gt.0) dodiis = .true.
!!!!===================================================================
!
!
!===================================================================
!         NEW VERSION
          if ( ndiis.gt.0 )  dodiis = .true.
!===================================================================
!
!
          if (dodiis) then
            xdiis(:,:,nmat) = sigma
            ediis(:,:,nmat) = sigma - sigold
            call diis(nbasis*nsph,nmat,xdiis,ediis,bmat,sigma)
          end if
!
!
!===================================================================
!         NEW VERSION
          do isph = 1,nsph
!
!           compute error
            delta(:) = sigma(:,isph) - sigold(:,isph)
!
!           energy norm of error      
            call hsnorm( delta(:), norm(isph) )
!
          end do

          drmsold = drms
          call rmsvec( nsph, norm, drms, dmax )
!!!          if ( drms > drmsold)  write(iout,1030)
!===================================================================
!
!
!         printing
          if (iprint.gt.1) then
            write(iout,1000) it, zero, drms, dmax
          end if
!
!         tolerance has been reached, exit loop
          if (drms.le.tol) goto 910
!
!         tolerance has NOT been reached, update solution
          sigold = sigma
!          
!       end of Jacobi iterations
        end do
!
!       UNSUCCESSFUL solve, end here
        write(iout,1020) 
        stop
!
!       SUCCESSFUL solve, just continue
  910   continue
!
      end if
      !
      ! free the memory:
      !
      if (iprint.gt.1) write(iout,*)
      if (star) deallocate(xi)
      deallocate( g, pot, vlm, sigold , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [1] failed deallocation !'
        stop
      endif
    !
      deallocate( delta, norm , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [2] failed deallocation !'
        stop
      endif
    !
      deallocate( vplm, basloc, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [3] failed deallocation !'
        stop
      endif
    !
      deallocate( xdiis, ediis, bmat , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv : [4] failed deallocation !'
        stop
      endif
    !
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
!      
 1000 format(' energy at iteration ',i4,': ',f14.7,' error (rms,max):',2f14.7)
 1010 format(' the solution to the ddCOSMO ',a,'equations took ',f8.3,' seconds.')
 1020 format(' ddCOSMO did not converge! Aborting...')
! 
      return
! 
end subroutine old_itsolv
!
!---------------------------------------------------------------------------------
! Purpose : computes a line of the (sparse) matrix/vector product for ddCOSMO.
!
!                   n                            n-1
!     pot = S \sigma  = sum \omega   * S   \sigma    + g
!                   i    j        ij    ij       j
!
!---------------------------------------------------------------------------------
subroutine old_calcv( first, isph, g, pot, sigma, basloc, vplm, vcos, vsin )
!
      implicit none
      logical,                        intent(in)    :: first
      integer,                        intent(in)    :: isph
      real*8, dimension(ngrid),       intent(in)    :: g
      real*8, dimension(nbasis,nsph), intent(in)    :: sigma
      real*8, dimension(ngrid),       intent(inout) :: pot
      real*8, dimension(nbasis),      intent(inout) :: basloc, vplm
      real*8, dimension(lmax+1),      intent(inout) :: vcos, vsin
!
      integer :: ig, ij, jsph
      real*8  :: vij(3), sij(3)
      real*8  :: vvij, tij, xij, oij
!---------------------------------------------------------------------------------
!
!     initialize
      pot = g
!
!     just return pot = g when n=1
      if ( first ) return
!
!     loop over integration points of i-sphere
      do ig = 1, ngrid
!
!       integration point also belongs to a neighboring sphere
        if ( ui(ig,isph).lt.one ) then
!
!         loop over neighboring spheres
          do ij = inl(isph), inl(isph+1) - 1
!
!           j-sphere is the neighbor
            jsph = nl(ij)
!
!           build omega_ij
            vij  = csph(:,isph) + rsph(isph)*grid(:,ig) - csph(:,jsph)
            vvij = sqrt(dot_product(vij,vij))
            tij  = vvij/rsph(jsph) 
            sij  = vij/vvij
            xij  = fsw(tij,se,eta)
            if (fi(ig,isph).gt.one) then
              oij = xij/fi(ig,isph)
            else
              oij = xij
            end if
!
!           point vij is inside j-sphere
!           ----------------------------
            if ( tij.lt.one ) then
!
!           compute spherical harmonics at integration point 
            call ylmbas( sij, basloc, vplm, vcos, vsin )
!                    
              pot(ig) = pot(ig) + oij*intmlp(tij,sigma(:,jsph),basloc)
!
!           point is on boundary of j-sphere (+ transition layer) [EXTENSION]
!           -----------------------------------------------------------------
            elseif ( tij.lt.( one + (se+one)/two*eta ) ) then
!
!             extension of potential
              select case(ext0)
!
!             t^l extension
              case(0)
!
!             compute spherical harmonics at integration point 
              call ylmbas( sij, basloc, vplm, vcos, vsin )
              pot(ig) = pot(ig) + oij*intmlp(tij,sigma(:,jsph),basloc)
!
!             constant extension
              case(1)
              pot(ig) = pot(ig) + oij*cstmlp(    sigma(:,jsph),basloc)
!
              endselect
!
            endif
!            
          end do
        end if
      end do
!
      return
!
end subroutine old_calcv
!
subroutine old_itsolv2( star, iefpcm, phi, psi, sigma, ene )
!
      implicit none
!      
      logical,                        intent(in)    :: star, iefpcm
      real*8, dimension(nbasis,nsph), intent(in)    :: phi, psi
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
      integer :: it, isph, nmat, lenb, ig, c1, c2, cr, istatus
      real*8  :: tol, drms, dmax, fep
      logical :: dodiis, first
      !
      integer, parameter :: nitmax=300
      real*8,  parameter :: ten=10.d0, tredis=1.0d-2
      !
      real*8 :: Lsigma(nbasis), s1, s2, s3, drmsold
      integer :: j
!      
!---------------------------------------------------------------------
!
!     initialize the timer:
      call system_clock(count_rate=cr)
      call system_clock(count=c1)
!
!     allocate local variables and set convergence parameters
      fep = one
      if (.not. iefpcm) fep = (eps-one)/eps
      tol = ten**(-iconv)
      allocate (g(ngrid,nsph),pot(ngrid),vlm(nbasis),sigold(nbasis,nsph) , &
      stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: allocation failed! [1]'
        stop
      endif

      sigold = zero
      allocate (delta(nbasis),norm(nsph) , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: allocation failed! [2]'
        stop
      endif

      allocate(vplm(nbasis),basloc(nbasis),vcos(lmax+1),vsin(lmax+1) , &
      stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: allocation failed! [3]'
        stop
      endif

      istatus=0
      if (star) allocate(xi(ngrid,nsph) , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: allocation failed! [4]'
        stop
      endif

      memuse = memuse + ngrid*nsph + ngrid*nproc + 2*nbasis*nproc + nbasis*nsph + &
               2*nbasis*nproc + 2*(lmax+1)*nproc + nsph
      if (star) memuse = memuse + nsph*ngrid
      memmax = max(memmax,memuse)
!
!     set up diis:
      dodiis = .false.
      nmat   = 1
      lenb   = ndiis + 1
      allocate (xdiis(nbasis,nsph,ndiis),ediis(nbasis,nsph,ndiis),bmat(lenb*lenb) , &
      stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: allocation failed! [5]'
        stop
      endif

      memuse = memuse + 2*nbasis*nsph*ndiis + lenb*lenb
      memmax = max(memmax,memuse)
!
!     solve the direct equations
!
      if (.not. star) then
!              
!       print header
        if (iprint.gt.1) then
          write(iout,*)'Solution of DIRECT system:'      
        endif
!
!
!!!!       build g:
!!!        g(:,:) = zero
!!!        call wghpot( phi, g )
!        
!       Jacobi iterations
        do it = 1, nitmax
!        
          first = ( it.eq.1 )
          vlm   = zero
!          
      !$omp parallel do default(shared) private(basloc,vcos,vsin,vplm) &
      !$omp private(isph,pot,vlm,delta) schedule(dynamic,10) 
!      
          do isph = 1, nsph
!          
            call calcv2( first, isph, pot, sigold, basloc, vplm, vcos, vsin )
            call intrhs( isph, pot, vlm )
!            
            vlm(:) = vlm(:) + phi(:,isph)
!            
            call solve( isph, vlm, sigma(:,isph) )
!            
!            
!!!!===================================================================
!!!            OLD VERSION
!!!            delta = sigma(:,isph) - sigold(:,isph)
!!!            call hsnorm(delta,norm(isph))
!!!!===================================================================
!            
!            
          enddo
!          
      !$omp end parallel do
!      
!            
!!!!===================================================================
!!!!         OLD VERSION
!!!          call rmsvec(nsph,norm,drms,dmax)
!!!          if ( drms.le.tredis .and. ndiis.gt.0 )  dodiis = .true.
!!!!===================================================================
!
!
!===================================================================
!         NEW VERSION
          if ( ndiis.gt.0 )  dodiis = .true.
!===================================================================
!
!          
          if (dodiis) then
!                  
            xdiis(:,:,nmat) = sigma
            ediis(:,:,nmat) = sigma - sigold
!            
            call diis(nbasis*nsph,nmat,xdiis,ediis,bmat,sigma)
!            
          end if
!          
!          
!===================================================================
!         NEW VERSION
          do isph = 1,nsph
!
!           compute error
            delta(:) = sigma(:,isph) - sigold(:,isph)
!
!           energy norm of error      
            call hsnorm( delta(:), norm(isph) )
!
          end do

          drmsold = drms
          call rmsvec( nsph, norm, drms, dmax )
!===================================================================
!
!
!         printing
          if ( iprint.gt.1 ) then
            ene = pt5*fep*sprod(nbasis*nsph,sigma,psi)
            write(iout,1000) it, ene, drms, dmax
          end if
!          
!         tolerance has been reached, exit loop
          if ( drms.le.tol ) goto 900
!
!         tolerance has NOT been reached, update solution
          sigold = sigma
!          
!       end of Jacobi iterations
        end do
!
!       SUCCESSFUL solve, just continue
        write(iout,1020) 
        stop
!        
!       SUCCESSFUL solve, just continue
  900   continue
!  
!
!       check solution
!       ==============
        if (iprint.gt.1) then
!
          s1 = zero ; s2 = zero ; s3 = zero ; pot = zero ; Lsigma = zero
!          
          do isph = 1,nsph
!
!           action of L_ij , j \ne i
            call calcv2( first, isph, pot, sigma, basloc, vplm, vcos, vsin )
            call intrhs( isph, pot, Lsigma )
            Lsigma = -Lsigma
!
!           action of L_ii
            Lsigma = Lsigma + sigma(:,isph) / facl
!            
            do j = 1,nbasis
!
              s1 = s1 + ( phi(  j,isph) - Lsigma(j) )**2
              s2 = s2 + ( phi(  j,isph)             )**2
              s3 = s3 + ( sigma(j,isph)             )**2
!              
            enddo
          enddo
!
          write(*,*) 'ddPCM :'
          write(*,*) ' '
          write(*,1001) sqrt(s3) , sqrt(s1) / sqrt(s2)
 1001     format(' || sigma || , || Phi_eps - L sigma || / || Phi_eps || = ',2(e12.5,2x) )       
!
        endif
!
!
!       compute energy : 1/2 f(eps) < psi , sigma >
        ene = pt5*fep*sprod(nbasis*nsph,sigma,psi)
!
      else
!              
!       print header
        if (iprint.gt.1) then
          write(iout,*)'Solution of ADJOINT system:'     
        endif
!
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
!
!
!!!!===================================================================
!!!!           OLD VERSION
!!!            delta = sigma(:,isph) - sigold(:,isph)
!!!            call hsnorm(delta,norm(isph))
!!!!===================================================================
!
!
          end do
      !$omp end parallel do
!
!
!!!!===================================================================
!!!!         OLD VERSION
!!!          call rmsvec(nsph,norm,drms,dmax)
!!!          if (drms.le.tredis .and. ndiis.gt.0) dodiis = .true.
!!!!===================================================================
!
!
!===================================================================
!         NEW VERSION
          if ( ndiis.gt.0 )  dodiis = .true.
!===================================================================
!
!
          if (dodiis) then
            xdiis(:,:,nmat) = sigma
            ediis(:,:,nmat) = sigma - sigold
            call diis(nbasis*nsph,nmat,xdiis,ediis,bmat,sigma)
          end if
!
!
!===================================================================
!         NEW VERSION
          do isph = 1,nsph
!
!           compute error
            delta(:) = sigma(:,isph) - sigold(:,isph)
!
!           energy norm of error      
            call hsnorm( delta(:), norm(isph) )
!
          end do

          drmsold = drms
          call rmsvec( nsph, norm, drms, dmax )
!!!          if ( drms > drmsold)  write(iout,1030)
!===================================================================
!
!
!         printing
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

      istatus=0
      if (star) deallocate(xi , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: deallocation failed! [1]'
        stop
      endif

      deallocate (g,pot,vlm,sigold , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: deallocation failed! [2]'
        stop
      endif

      deallocate (delta,norm , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: deallocation failed! [3]'
        stop
      endif

      deallocate (vplm,basloc,vcos,vsin , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: deallocation failed! [4]'
        stop
      endif

      deallocate (xdiis,ediis,bmat , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'itsolv2: deallocation failed! [5]'
        stop
      endif

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
 1000 format(' energy at iteration ',i4,': ',f14.7,' error(rms,max):',2(e12.5,2x))
 1010 format(' the solution to the ddCOSMO ',a,'equations took ',f8.3,' seconds.')
 1020 format(' ddCOSMO did not converge! Aborting...')
      return
!
end subroutine old_itsolv2
!
