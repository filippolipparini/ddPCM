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
!------------------------------------------------------------------------------      
! Compute contraction :
!
!   f =+ < s , dA/dr_i x > = 
!
!     =+ sum s_j ( sum (d_i A)_jk x_k )
!         j         k
!
! Recall that (d_i A)_jj = 0 for i \ne j . Then :
!
!     =+ sum s_j (   sum  (d_i A)_jk x_k + (d_i A)_jj x_j )
!         j        k \ne j
!
!     =+ sum s_j    sum  (d_i A)_jk x_k + sum s_j (d_i A)_jj x_j
!         j       k \ne j                  j
!
!     =+ sum s_j    sum  (d_i A)_jk x_k + s_i (d_i A)_ii x_i
!         j       k \ne j 
!
!     =+   sum   s_j   sum   (d_i A)_jk x_k + s_i   sum   (d_i A)_ik x_k +  s_i (d_i A)_ii x_i
!        j \ne i     k \ne j                      k \ne i
!
! Finally, since (d_i A)_jk = 0 for k \ne i and j \ne i , then :
!
!     =+   sum   s_j (d_i A)_ji x_i + s_i   sum   (d_i A)_ik x_k +  s_i (d_i A)_ii x_i
!        j \ne i                          k \ne i
!
!             +---+                       +--------------------+        +---+
!             |///|                       |///////////|  |/////|        |///|
!             |///|                       +--------------------+        +---+
!             |///|
!             |///|
!             |///|
!             |///|
!             |---|
!             |   |
!             |---|
!             |///|
!             +---+
!------------------------------------------------------------------------------      
subroutine service_routine1( s, x, isph, f )
!
      use ddcosmo , only : ui, nsph, nbasis, zero, ngrid, w, one, basis, csph, &
                           rsph, grid, zi, lmax, dbasis, du
!
      implicit none 
      real*8, dimension(nbasis,nsph), intent(in)    :: s
      real*8, dimension(nbasis,nsph), intent(in)    :: x
      integer,                        intent(in)    :: isph
      real*8, dimension(3),           intent(inout) :: f
!
      real*8,  dimension(nbasis) :: basloc, vplm
      real*8,  dimension(3,nbasis) :: dbsloc
      real*8,  dimension(lmax+1) :: vcos, vsin
      real*8, dimension(3) :: f4,vij,s_ijn,dt_ijn,f2
      real*8, dimension(3,3) :: ds_ijn
      real*8 :: vvij,t_ijn,f3!!,err
      integer :: n,icomp,jcomp,jsph,ksph
!      
!------------------------------------------------------------------------------      
!
!     case j \ne i : I_1
!     ==================
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize f4(n)
        f4(1:3) = zero
!
!       1st loop over spheres
        do jsph = 1,nsph
!
          if ( jsph.ne.isph ) then
!
!           non-zero contribution
            if ( ui(n,jsph).gt.zero ) then
!
!             compute f3(j,n)
!             ---------------
              f3 = ui(n,jsph) * dot_product( basis(:,n), s(:,jsph) )
!
!             compute s_ijn, t_ijn
              vij(:)   = csph(:,jsph) + rsph(jsph)*grid(:,n) - csph(:,isph)
              vvij     = sqrt( dot_product( vij(:), vij(:) ) )
              t_ijn    = rsph(isph)/vvij 
              s_ijn(:) =     vij(:)/vvij
!
!             compute derivatives \grad_i t_ijn
              dt_ijn(1:3) = rsph(isph) * vij(1:3) / vvij**3
!
!             compute derivatives ds_ijn,icomp / dr_i,jcomp = ds_ijn(icomp,jcomp)
              do icomp = 1,3
                do jcomp = 1,3
!
                  ds_ijn(icomp,jcomp) = vij(icomp)*vij(jcomp) / vvij**3
!
                  if ( icomp.eq.jcomp ) then
!
                    ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) - one / vvij
!
                  endif
!                  
                enddo
              enddo

!!!              err=zero
!!!              do icomp=1,3
!!!              do jcomp=1,3
!!!                err = err + (ds_ijn(icomp,jcomp)-ds_ijn(jcomp,icomp))**2
!!!              enddo
!!!              enddo
!!!              write(*,*)'err = ',sqrt(err)

!             
!             compute Y_ll^mm(s_ijn)
              call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!             compute f2(j,n)
!             ---------------
              call compute_f2_jn( t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,isph), f2(1:3) )
!                      
!             accumulate for f4(n)
!             --------------------
              f4(1:3) = f4(1:3) + f2(1:3) * f3
!                      
            endif  
!
!
!           missing term
!           ============
!
!           non-zero contribution from derivative d_i ( U_j^n )
            if ( sum(abs( du(:,isph,n,jsph) )).gt.zero ) then
!
!             d_i U_j^n * Y_l^m(s_n) * s_j          
              f2(1:3) = du(1:3,isph,n,jsph) * dot_product( basis(:,n), s(:,jsph) )
!
!             compute f2(n) [ store in f3 for convenience ]
!             ---------------------------------------------
              call compute_f2_n( basis(:,n), x(:,isph), f3 )
!
!             accumulate for f4
!             ----------------
              f4(1:3) = f4(1:3) - f2(1:3)*f3
!!!              write(*,*)'cont = ',f2*f3
!!!              read(*,*)
!
            endif
!
!
          endif 
!
        enddo  
!
!       accumulate for f
!       ----------------
        f(1:3) = f(1:3) - w(n) * f4(1:3)
!
      enddo  
!
!
!     case j \eq i : I_2,1
!     ====================
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize f4(n)
        f4(1:3) = zero
!
!       1st loop over spheres
        do ksph = 1,nsph
!
          if ( ksph.ne.isph ) then
!
!           compute f3(k,n)
!           ---------------
            f3 = dot_product( basis(:,n), s(:,isph) )
!
!           compute s_kjn, t_kjn
            vij    = csph(:,isph) + rsph(isph)*grid(:,n) - csph(:,ksph)
            vvij   = sqrt(dot_product(vij,vij))
            t_ijn  = rsph(ksph)/vvij 
            s_ijn  =        vij/vvij
!
!           compute derivatives \grad_j t_kjn
            dt_ijn(1:3) = - rsph(ksph) * vij(1:3) / vvij**3
!
!           compute derivatives ds_kjn,icomp / dr_j,jcomp = ds_ijn(icomp,jcomp)
            do icomp = 1,3
              do jcomp = 1,3
!
                ds_ijn(icomp,jcomp) = - vij(icomp)*vij(jcomp) / vvij**3
!
                if ( icomp.eq.jcomp ) then
!
                  ds_ijn(icomp,jcomp) = ds_ijn(icomp,jcomp) + one / vvij
!
                endif
!                
              enddo
            enddo
!           
!           compute Y_ll^mm(s_kjn)
            call dbasis( s_ijn, basloc, dbsloc, vplm, vcos, vsin )
!
!           compute f2(j,n)
!           ---------------
            call compute_f2_kn( isph, n, t_ijn, dt_ijn, ds_ijn, basloc, dbsloc, x(:,ksph), f2(1:3) )
!                    
!           accumulate for f4(n)
!           --------------------
            f4(1:3) = f4(1:3) + f2(1:3) * f3
!                      
          endif  
        enddo  
!
!       accumulate for f
!       ----------------
        f(1:3) = f(1:3) - w(n) * f4(1:3)
!
      enddo  
!
!
!     case j \eq i : I_2,2
!     ====================
!
!     loop over integration points
      do n = 1,ngrid
!
!       compute f3(n) [ store in f2 for convenience ]
!       ---------------------------------------------
        f2(1:3) = zi(1:3,n,isph) * dot_product( basis(:,n), s(:,isph) )
!
!       compute f2(n) [ store in f3 for convenience ]
!       ---------------------------------------------
!!!        call compute_f2_n( basloc, x(:,isph), f3 )
!!!        write(*,*)'old f3 = ',f3
        call compute_f2_n( basis(:,n), x(:,isph), f3 )
!!!        write(*,*)'new f3 = ',f3
!!!        read(*,*)
!
!       accumulate for f
!       ----------------
        f(1:3) = f(1:3) + w(n) * f2(1:3) * f3
!
      enddo  
!
!
endsubroutine service_routine1
!----------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------
subroutine compute_f2_jn( t, dt, ds, basloc, dbsloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, four, pi
!
      implicit none
      real*8,                      intent(in)  :: t
      real*8, dimension(3  ),      intent(in)  :: dt
      real*8, dimension(3,3),      intent(in)  :: ds
      real*8, dimension(  nbasis), intent(in)  :: basloc
      real*8, dimension(3,nbasis), intent(in)  :: dbsloc
      real*8, dimension(  nbasis), intent(in)  :: x
      real*8, dimension(3),        intent(out) :: f2
!
      real*8, dimension(3,nbasis) :: s1,s2,f1
      real*8, dimension(3) :: s3
      real*8 :: tt,fl,fac
      integer :: icomp,jcomp,l,m,ind
!      
!----------------------------------------------------------------------------------
! Recall :
!
!   f1(1:3,lm) = \grad [ t^l+1 * Y(s) ]
!
!              = t^l  [ (l+1) \grad t Y(s) + t ( \grad s )^T \grad Y(s) ]
!
!              = tt   [ (l+1) s1(1:3,lm)   + s2(1:3,lm)                 ]
!
! and :
!
!                  4 pi l
!   f2(1:3) = sum  ------  sum x(lm) * f1(1:3,lm)
!              l   2l + 1   m
!
!                  4 pi l
!           = sum  ------  s3(1:3,l)  
!              l   2l + 1 
!
!           = sum fac(l) * s3(1:3,l) 
!              l
!
!----------------------------------------------------------------------------------
!
!     initialize f2
      f2(1:3) = zero
!
!     compute s1
      do icomp = 1,3
!      
        s1(icomp,1:nbasis) = dt(icomp)*basloc(1:nbasis)
!        
      enddo
!
!     compute s2
      do icomp = 1,3
!
        s2(icomp,1:nbasis) = zero
!
!       accumulate
        do jcomp = 1,3
!
          s2(icomp,1:nbasis) = s2(icomp,1:nbasis) + dbsloc(jcomp,1:nbasis)*ds(jcomp,icomp)
!
        enddo
!
        s2(icomp,1:nbasis) = t * s2(icomp,1:nbasis) 
!
      enddo
!
!     initialize factor t^l
      tt = one
!
!     contract over l
      do l = 0,lmax
!
!       build factor : 4pi*l / (2l+1)
        fl = dble(l)
        fac = four*pi*fl / (two*fl+one)
!
!       build f1
        f1(1:3,1:nbasis) = tt * ( (fl+one)*s1(1:3,1:nbasis) + s2(1:3,1:nbasis) )
!
!       compute 1st index
        ind = l*l + l + 1
!
!       initialize s3(l)
        s3(1:3) = zero
!
!       contract over m
        do m = -l,l
!
          s3(1:3) = s3(1:3) + f1(1:3,ind+m) * x(ind+m)          
!
        enddo
!
!       accumulate for f2
        f2(1:3) = f2(1:3) + fac*s3(1:3)
!
!       update factor t^l
        tt = tt*t
!
      enddo
!
!
endsubroutine compute_f2_jn
!----------------------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------------------
subroutine compute_f2_kn( jsph, n, t, dt, ds, basloc, dbsloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, four, pi, ui, zi
!
      implicit none
      integer,                     intent(in)  :: jsph
      integer,                     intent(in)  :: n
      real*8,                      intent(in)  :: t
      real*8, dimension(3  ),      intent(in)  :: dt
      real*8, dimension(3,3),      intent(in)  :: ds
      real*8, dimension(  nbasis), intent(in)  :: basloc
      real*8, dimension(3,nbasis), intent(in)  :: dbsloc
      real*8, dimension(  nbasis), intent(in)  :: x
      real*8, dimension(3),        intent(out) :: f2
!
      real*8, dimension(3,nbasis) :: s1,s2,f1
      real*8, dimension(3) :: s3
      real*8 :: tt,fl,fac
      integer :: icomp,jcomp,l,m,ind
!      
!----------------------------------------------------------------------------------
! Recall :
!
!   f1(1:3,lm) = \grad [ U * t^l+1 * Y(s) ]
!
!              = (\grad U) t^l+1 * Y(s) + U * t^l [ (l+1) (\grad t) Y(s) + t ( \grad s )^T \grad Y(s) ]
!
!              = t^l [ (l+1) (\grad t) U * Y + t [ (\grad U) Y + U ( \grad s )^T \grad Y ] ]
!
!              = tt  [ (l+1) s1(1:3,lm)      + t   s2(1:3,lm) ]
!
! and :
!
!                  4 pi l
!   f2(1:3) = sum  ------  sum x(lm) * f1(1:3,lm)
!              l   2l + 1   m
!
!                  4 pi l
!           = sum  ------  s3(1:3,l)  
!              l   2l + 1 
!
!           = sum fac(l) * s3(1:3,l) 
!              l
!
!----------------------------------------------------------------------------------
!
!     initialize f2
      f2(1:3) = zero
!
!     compute s1
      do icomp = 1,3
!      
        s1(icomp,1:nbasis) = dt(icomp) * basloc(1:nbasis) * ui(n,jsph)
!        
      enddo
!
!     compute s2
      s2(1:3,1:nbasis) = zero
!      
      do icomp = 1,3
!
!       accumulate
        do jcomp = 1,3
!
          s2(icomp,1:nbasis) = s2(icomp,1:nbasis) + dbsloc(jcomp,1:nbasis)*ds(jcomp,icomp)
!
        enddo
!
        s2(icomp,1:nbasis) = ui(n,jsph)*s2(icomp,1:nbasis) + zi(icomp,n,jsph)*basloc(1:nbasis)
!
      enddo
!
!     initialize factor t^l
      tt = one
!
!     contract over l
      do l = 0,lmax
!
!       build factor : 4pi*l / (2l+1)
        fl = dble(l)
        fac = four*pi*fl / (two*fl+one)
!
!       build f1
        f1(1:3,1:nbasis) = tt * ( (fl+one)*s1(1:3,1:nbasis) + t*s2(1:3,1:nbasis) )
!
!       compute 1st index
        ind = l*l + l + 1
!
!       initialize s3(l)
        s3(1:3) = zero
!
!       contract over m
        do m = -l,l
!
          s3(1:3) = s3(1:3) + f1(1:3,ind+m) * x(ind+m)          
!
        enddo
!
!       accumulate for f2
        f2(1:3) = f2(1:3) + fac*s3(1:3)
!
!       update factor t^l
        tt = tt*t
!
      enddo
!
!
endsubroutine compute_f2_kn
!----------------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------------
subroutine compute_f2_n( basloc, x, f2 )
!
      use ddcosmo , only : lmax, nbasis, zero, one, two, pi
!
      implicit none
      real*8, dimension(nbasis), intent(in)  :: basloc
      real*8, dimension(nbasis), intent(in)  :: x
      real*8,                    intent(out) :: f2
!
      real*8 :: fl,s,fac
      integer :: l,m,ind
!      
!----------------------------------------------------------------------------------
! Recall :
!
!             2 pi 
!   f2 = sum  ----  sum Y(lm) * x(lm)
!         l   2l+1   m
!
!             2 pi
!      = sum  ----  s(l)
!         l   2l+1 
!
!----------------------------------------------------------------------------------
!
!     initialize f2
      f2 = zero
!
!     contract over l
      do l = 0,lmax
!
!       build factor : 2pi / (2l+1)
        fl = dble(l)
        fac = two*pi/(two*fl+one)
!
!       compute 1st index
        ind = l*l + l + 1
!
!       initialize s
        s = zero
!        
!       contract over m
        do m = -l,l
!
          s = s + basloc(ind+m) * x(ind+m)          
!
        enddo
!
!       accumulate for f2
        f2 = f2 + fac*s
!
      enddo
!
!
endsubroutine compute_f2_n
!----------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
real*8 function extmlp(t,sigma,basloc)
implicit none
real*8, intent(in) :: t
real*8, dimension(nbasis), intent(in) :: sigma, basloc
!
integer :: l, ind
real*8  :: tt, ss, fac
!
tt = one/t
ss = zero
do l = 0, lmax
  ind = l*l + l + 1
  fac = tt/facl(ind)
  ss = ss + fac*dot_product(basloc(ind-l:ind+l),sigma(ind-l:ind+l))
  tt = tt/t
end do
extmlp = ss
return
endfunction extmlp
!---------------------------------------------------------------------
!
!
!
!!!!------------------------------------------------------------------------------------------------
!!!real*8 function fsw( t, s, eta )
!!!!
!!!! Remark : a step-function symmetrically squeeshing around x = 1
!!!!
!!!! The "s" variable is totally redundant. In fact :
!!!!
!!!!   x = t - eta*s  ==> y = x - 0.5*eta = t - (0.5 + s)*eta
!!!!
!!!! Thus :
!!!!
!!!!   x >= f_hi  ==>  t >= 1 + (0.5 + s)*eta
!!!!   . . .           t <= 1 - (0.5 + s)*eta
!!!!
!!!
!!!!      
!!!!  s = \hat{s} eta ; eta \in [0,1]
!!!!     
!!!
!!!!------------------------------------------------------------------------------------------------
!!!!
!!!      implicit none
!!!      real*8, intent(in) :: t, s, eta
!!!!      
!!!      real*8 :: a, b, x, y, f_hi, f_low
!!!      real*8, parameter :: zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0, f6=6.0d0, f10=10.d0, &
!!!                           f12=12.d0, f15=15.d0
!!!!                           
!!!!------------------------------------------------------------------------------------------------
!!!!
!!!      x = t - eta*s 
!!!!
!!!!     auxiliary variable for function \chi [ not really needed !!! ]
!!!      y = x - pt5*eta
!!!!
!!!!     lower and upper bounds of transition area of function \chi
!!!      f_low = one - pt5*eta
!!!      f_hi  = one + pt5*eta
!!!!      
!!!!     construct smooth step function \chi(x)
!!!      if     ( x.ge.f_hi  ) then
!!!!              
!!!        fsw = zero
!!!!        
!!!      elseif ( x.le.f_low ) then
!!!!      
!!!        fsw = one
!!!!        
!!!      else
!!!!              
!!!        fsw = ( (y-one)*(y-one) * (y-one+two*eta)*(y-one+two*eta) ) / (eta**4)
!!!!              
!!!      endif
!!!!      
!!!      return
!!!!
!!!!
!!!endfunction fsw
!!!!------------------------------------------------------------------------------------------------
!
!
!!!!------------------------------------------------------------------------------------------------
!!!!     switching function derivative for ddCOSMO regularization.
!!!!
!!!real*8 function dfsw( t, s, eta )
!!!!
!!!      implicit none
!!!      real*8, intent(in) :: t, eta, s
!!!!
!!!      real*8  flow, fhi, x, y
!!!      real*8, parameter :: f30=30.0d0
!!!!      
!!!!------------------------------------------------------------------------------------------------
!!!!
!!!      x = t - eta*s 
!!!      y = x - pt5*eta
!!!!
!!!!!!      flow = one - eta
!!!      flow = one - pt5*eta
!!!      fhi  = one + pt5*eta
!!!!      
!!!!!!      if     ( t.ge.one ) then
!!!      if     ( x.ge.fhi ) then
!!!!              
!!!        dfsw = zero
!!!!        
!!!      elseif ( x.le.flow ) then
!!!!        
!!!!!!        dfsw = one
!!!        dfsw = zero
!!!!        
!!!      else
!!!!        
!!!!!!        dfsw = f30*(one-t)*(t-one)*(t-one+eta)*(t-one+eta) / (eta**5)
!!!        dfsw = 4.d0 * (y - 1.d0 - 0.5d0*eta) * &
!!!                      (y - 1.d0 + 0.5d0*eta) * &
!!!                      (y - 1.d0 + 1.5d0*eta) / eta**4
!!!!        
!!!      endif
!!!!        
!!!      return
!!!!        
!!!!        
!!!endfunction dfsw
 
