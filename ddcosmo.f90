!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
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
! grids by D. Laikov and C. van Wuellen, as publicly available on CCL. If the
! routines are used, the following reference should also be included:
!
! [6] V.I. Lebedev, and D.N. Laikov
!     "A quadrature formula for the sphere of the 131st
!      algebraic order of accuracy"
!     Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
! Written by Filippo Lipparini, October 2015.
!
!-------------------------------------------------------------------------------
!
module ddcosmo
!
      implicit none
!
!     - arguments contained in bash script
!
      character(len=64), dimension(2) :: args
!
!     - numerical constants
!
      integer, parameter :: ndiis=25, iout=6, nngmax=100
      real*8,  parameter :: zero=0.d0, pt5=0.5d0, one=1.d0, two=2.d0, four=4.d0
      real*8,  parameter :: toang=0.52917721092d0, tokcal=627.509469d0
!
!     - numerical constants explicitly computed
!
      real*8  :: pi, sq2
!
!     - quantities contained in control file
!
!     iprint     - printing flag
!     nproc      - number of openMP threads ; 0) no parallelization
!     lmax       - max angular momentum of spherical harmonics basis
!     ngrid      - desired number of Lebedev integration points
!     iconv      - threshold for iterative solver ( 10^-iconv )
!     igrad      - 1) compute forces ; 0) do not compute forces
!     iscrf      - 0) use cosmo ; 1) use pcm
!     eps        - dielectric constant of the solvent
!     iunit      - 0) convert to bohr ; 1) do not convert to bohr
!     eta, se    - regularization parameters
!     ext0, ext1 - extension of potential for COSMO and PCM
!     isolver    - 0) jacobi/diis ; 1) gmres
!
      integer :: iprint, nproc, lmax, ngrid, iconv, igrad, iscrf, iunit, ext0, ext1, isolver
      real*8  :: eps, eta, se           
!
!     - other quantities
!
!     nsph   - number of spheres/atoms
!     ncav   - number of integration points on cavity's boundary
!     nbasis - number of basis functions, i.e., spherical harmonics
!
      integer :: nsph, ncav, nbasis
!
!     - memory usage
!
      integer :: memuse, memmax
!
!     - workspaces
!
      integer, allocatable :: inl(:), nl(:)
      real*8,  allocatable :: rsph(:), csph(:,:), ccav(:,:)
      real*8,  allocatable :: w(:), grid(:,:), basis(:,:)
      real*8,  allocatable :: fact(:), facl(:), facs(:)
      real*8,  allocatable :: fi(:,:), ui(:,:), zi(:,:,:), du(:,:,:,:)
      real*8,  allocatable :: prec(:,:,:), precm1(:,:,:)
      real*8,  allocatable :: read_x(:),read_y(:),read_z(:),read_r(:),read_q(:)
!
!     - Lebedev integration rules
!    
!     ngrid_vec(i) - number of integration points in i-th grid
!     lmax_vec( i) - angular momentum integrated by i-th grid
!     nLLG         - number of grids supported by the code      
!
      integer, parameter, dimension(32) :: ngrid_vec = (/   6,  14,  26,  38,  50,  74,  86, 110,  &
                                                          146, 170, 194, 230, 266, 302, 350, 434,  &
                                                          590, 770, 974,1202,1454,1730,2030,2354,  &
                                                         2702,3074,3470,3890,4334,4802,5294,5810  /)
      integer, parameter, dimension(32) :: lmax_vec  = (/   3,   5,   7,   9,  11,  13,  15,  17,  &
                                                           19,  21,  23,  25,  27,  29,  31,  35,  &
                                                           41,  47,  53,  59,  65,  71,  77,  83,  &
                                                           89,  95, 101, 107, 113, 119, 125, 131  /)
      integer, parameter :: nLLG = 32
!
!     - miscellanea
!
      logical :: grad, do_diag
      logical :: iquiet = .false.
      logical :: use_fmm = .false.
!
!     - subroutines & functions
!
!     read_control_file  - read control file
!     read_molecule_file - read molecule file
!     reset_ngrid        - choose suitable grid of integration points
!     set_pi             - compute numerical constants
!     ddinit             - initialize data structure
!     memfree            - free data structure
!     sprod              - scalar product
!     fsw                - smoothing function
!     dfsw               - derivative of smoothing function
!     ptcart             - print routine
!     prtsph             - print routine
!     intrhs             - integrate spherical harmonics expansions
!     ylmbas             - compute spherical harmonics Y_l^m
!     dbasis             - compute derivatives of spherical harmonics
!     polleg             - compute Legendre polynomials
!     trgev              - service routine for computation of spherical harmonics
!     intmlp             -
!     wghpot             -
!     hsnorm             -
!     adjrhs1            - auxiliary routine for COSMO adjoint action
!     header             - print header
!     fdoka              -
!     fdokb              -
!     fdoga              -
!     calcv2             - auxiliary routine for COSMO action
!
      contains
!
!
!--------------------------------------------------------------------------------------------------
subroutine read_control_file()
!
      implicit none
!      
!--------------------------------------------------------------------------------------------------
!
!     open control file
      open( unit=10, file=args(1) )

!     read control parameters
      read(10,*) iprint      ! printing flag
      read(10,*) nproc       ! number of openmp threads
      read(10,*) lmax        ! max angular momentum of spherical harmonics basis
      read(10,*) ngrid       ! number of lebedev points
      read(10,*) iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
      read(10,*) igrad       ! whether to compute (1) or not (0) forces
      read(10,*) iscrf       ! whether to use cosmo (0) or pcm (1)
      read(10,*) eps         ! dielectric constant of the solvent
      read(10,*) iunit       ! whether to convert to bohr (0) or not (1)
      read(10,*) eta, se     ! regularization parameters
      read(10,*) ext0, ext1  ! extension of potential for COSMO and PCM
      read(10,*) isolver     ! whether to use the jacobi/diis (0) or gmres (1) solver
!
!     close control file
      close(10)
!
!
endsubroutine read_control_file
!--------------------------------------------------------------------------------------------------
!
!
!
!
!--------------------------------------------------------------------------------------------------
subroutine read_molecule_file()
!
      implicit none
!
      real*8 :: tobohr
      integer :: i,istatus
!      
!--------------------------------------------------------------------------------------------------
!
!     open molucule file
      open( unit=10, file=args(2) )
!
!     read number of atoms
      read(10,*) nsph
!
!     allocate arrays for centers, radii, charges
      allocate( read_x(nsph), read_y(nsph), read_z(nsph), read_r(nsph), read_q(nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'read_molecule_file : failed allocation !'
        stop
      endif
!
!     update memory usage
      memuse = memuse + 5*nsph
      memmax = max(memmax,memuse)
!
!     read charges, centers, radii
      do i = 1, nsph
        read(10,*) read_q(i), read_x(i), read_y(i), read_z(i), read_r(i)
      enddo
!      
!     if required, convert to Angstrom 
      if ( iunit.eq.0 ) then
!              
        tobohr = 1.0d0/toang
!      
        read_x = read_x*tobohr
        read_y = read_y*tobohr
        read_z = read_z*tobohr
        read_r = read_r*tobohr
!        
      endif
!
!     close molecule file
      close(10)
!
!
endsubroutine read_molecule_file
!--------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------
!
subroutine reset_ngrid( igrid )
! 
      implicit none
      integer, intent(out) :: igrid
!
      integer :: iflag, ig, idec
!
!-----------------------------------------------------------------------------------
!
!!!      do lmax=2,80
!
!     initialize control flag
      iflag = 0
!
!     initialize to largest grid
      igrid=nLLG

!     loop over grids
      do ig=1,nLLG
!
!       grid can integrate 2*lmax angular momentum 
        if ( lmax_vec(ig).ge.2*lmax ) then
!                
!         save grid number
          igrid=ig
!
!         update control flag, break out of loop
          iflag = 1 ; exit
!          
        endif
!        
      enddo
!      
!     adjust ngrid
      ngrid = ngrid_vec(igrid)
!
!     check control flag
      if ( iflag.eq.0 ) then 
!              
        write(*,1000) lmax      
 1000   format(' Integration grid for lmax = ',i3,' not available !')
        write(*,*)'Largest grid selected. Continue ? 1 - Yes ; 0 - No'
        read(*,*) idec
        if ( idec.ne.1 )  stop
!        
      endif
!
!!!        write(*,*)'lmax,igrid,ngrid = ',lmax,igrid,ngrid
!!!      enddo
!!!      stop
!
!
endsubroutine reset_ngrid
!-----------------------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------------------
subroutine set_pi
!
      implicit none
!      
      pi  = four*atan(one)
      sq2 = sqrt(two)
!
!
endsubroutine set_pi
!----------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------
! Purpose : allocate the various arrays needed for ddcosmo, assemble the cavity
!           and the various associated geometrical quantities.
!----------------------------------------------------------------------------------
subroutine ddinit( n, x, y, z, rvdw )
!
      implicit none
      integer,               intent(in) :: n
      real*8,  dimension(n), intent(in) :: x, y, z, rvdw
!
      integer :: isph, jsph, i, ii, lnl, l, ind, m, igrid, inear, jnear, &
      istatus, ig, ji
      real*8  :: fac, fl, ffl, fnorm, d2, r2, v(3), vv, t, xt, swthr
!
      real*8,  allocatable :: vcos(:), vsin(:), vplm(:)
      real*8 :: tji,vji(3),vvji,sji(3)
!
!-------------------------------------------------------------------------------
!
!     openMP parallelization
      if ( nproc.eq.0 )  nproc = 1
!    
!    $ call omp_set_num_threads(nproc)
!    
!     compute numerical constants
      call set_pi
!    
!     print header
      if ( .not.iquiet ) call header
!    
!     compute forces flag
      grad = ( igrad.ne.0 )
!      
!     number of basis functions
      nbasis = (lmax+1)*(lmax+1)
!    
!     allocate quantities
      allocate( rsph(nsph), csph(3,nsph), w(ngrid), grid(3,ngrid), basis(nbasis,ngrid), &
                inl(nsph+1), nl(nsph*nngmax), fi(ngrid,nsph), ui(ngrid,nsph), &
                du(3,nsph,ngrid,nsph), zi(3,ngrid,nsph), fact(max(2*lmax+1,2)), &
                facl(nbasis), facs(nbasis) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'ddinit : [1] allocation failed !'
        stop
      endif
!    
!     update memory usage
      memuse = memuse + 4*nsph + 4*ngrid + nbasis*ngrid + nsph+1 + nsph*nngmax + &
               2*ngrid*nsph + 2*lmax+1 + 2*nbasis
      if ( grad )  memuse = memuse + 3*ngrid*nsph
      memmax = max(memmax,memuse)
!    
!     compute factorials
      fact(1) = one
      fact(2) = one
!      
      do i = 3,2*lmax+1
!
        fact(i) = dble(i-1)*fact(i-1)
!        
      enddo
!    
!     compute factors for spherical harmonics
      do l = 0,lmax
!      
        ind = l*l + l + 1
!        
        fl  = (two*dble(l) + one)/(four*pi)
        ffl = sqrt(fl)
        facl( ind-l:ind+l) = fl
        facs(ind) = ffl
!        
        do m = 1,l
!        
          fnorm = sq2*ffl*sqrt(fact(l-m+1)/fact(l+m+1))
          if ( mod(m,2).eq.1 ) fnorm = -fnorm
!          
          facs(ind+m) = fnorm
          facs(ind-m) = fnorm
!          
        enddo
      enddo
!    
!     set the centers and radii of the spheres
      csph(1,:) = x
      csph(2,:) = y
      csph(3,:) = z
      rsph      = rvdw
!    
!     load integration points grid
      call llgrid( ngrid, w, grid )
!    
!     allocate workspaces
      allocate( vplm(nbasis), vcos(lmax+1), vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'ddinit : [2] allocation failed !'
        stop
      endif
!    
!     update memory usage
      memuse = memuse + nproc*(nbasis + 2*lmax + 2)
      memmax = max(memmax,memuse)
!    
      !$omp parallel do default(shared) private(i,vplm,vcos,vsin)
!      
!     loop over integration points
      do i = 1,ngrid
!
!       compute spherical hamonics at grid point
        call ylmbas( grid(:,i), basis(:,i), vplm(:), vcos(:), vsin(:) )
!        
      enddo
!      
      !$omp end parallel do
!      
!     deallocate workspaces
      deallocate( vplm, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'ddinit : [1] deallocation failed!'
        stop
      endif
!    
!     update memory usage
      memuse = memuse - nproc*(nbasis + 2*lmax + 2)
!    
!    
!     build neighbors list (CSR format)
!     =================================
!
!      \\  jsph |
!     isph  \\  |  1   2   3   4   5   6
!     -----------------------------------
!             1 |      x       x   x
!             2 |  x       x       x   x
!             3 |      x       x       x
!             4 |  x       x       x   x
!             5 |  x   x       x
!             6 |      x   x   x        
!
!
!      inl =  1,       4,          8,      11,         15,      18,21        pointer to 1st neighbor
!
!             |        |           |        |           |        |
!             v        V           V        V           V        V
!
!             1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|16|17|18|19|20
!
!      nl  =  2, 4, 5, 1, 3, 5, 6, 2, 4, 6, 1, 3, 5, 6, 1, 2, 4, 2, 3, 4     neighbors list
!
!    
!     index of nl  
      ii = 1
!
!     number of neighbors recored in nl so far
      lnl = 0
!    
!     loop over i-spheres
      do isph = 1,nsph
!    
!       pointer to 1st neighbor of i-sphere
        inl(isph) = lnl + 1
!    
!       loop over j-spheres
        do jsph = 1, nsph
!    
!         exclude i-sphere from neighbors of i-sphere
          if ( isph.ne.jsph ) then
!    
!           distance square b/w centers
            d2 = (csph(1,isph) - csph(1,jsph))**2 + &
                 (csph(2,isph) - csph(2,jsph))**2 + &
                 (csph(3,isph) - csph(3,jsph))**2
!            
!           sum square of radii, accounting for switch region
            r2 = ( rsph(isph)*( 1.d0+(se+1.d0)*eta/2.d0 ) + &
                   rsph(jsph)*( 1.d0+(se+1.d0)*eta/2.d0 )   )**2
!    
!           spheres intersect
            if ( d2.le.r2 ) then
!    
!             record neighbor
              nl(ii) = jsph
!    
!             advance index of nl vector
              ii  = ii + 1
!
!             increment number of neighbors recorded so far
              lnl = lnl + 1
!              
            endif
          endif
        enddo
      enddo
!      
!     last entry for consistency with CSR format
      inl(nsph+1) = lnl+1
!    
!    
!-----------------------------------------------------------------------
! Define :
!
!   N_i = list of neighbors of i-sphere [ excluding i-sphere ]
!
!            | r_i + rho_i s_n - r_j |
!   t_n^ij = -------------------------
!                      rho_j
!
!   fi(n,i) =    sum    \chi( t_n^ij )
!             j \in N_i 
! 
! Notice that the derivative of fi(n,i) wrt to r_k is (possibly) nonzero
! when, either k = i, or k \in N_j .
!
! Define :
!
!             [  1 - fi(n,i)   ,   when fi(n,i) <= 1
!   ui(n,i) = [
!             [  0             ,   otherwise
!
! REMARK : zi(n,i) is a close relative to d ui(n,i) / dr_i , but not
!          quite the same ...
!-----------------------------------------------------------------------
!    
!     build arrays fi, ui, zi
!     =======================
!
!     initialize
      fi(:,:)=zero ; ui(:,:)=zero ; if ( grad )  zi(:,:,:)=zero
!      
      !$omp parallel do default(shared) private(isph,i,ii,jsph,v,vv,t,xt,swthr,fac)
!    
!     loop over spheres
      do isph = 1,nsph
!      
!       loop over integration points
        do i = 1,ngrid
!    
!         loop over neighbors of i-sphere
          do ii = inl(isph),inl(isph+1)-1
!    
!           neighbor's number
            jsph = nl(ii)
!            
!           compute t_n^ij
            v(:) = csph(:,isph) + rsph(isph)*grid(:,i) - csph(:,jsph)
            vv   = sqrt(dot_product(v,v))
            t    = vv/rsph(jsph)
!    
!           compute \chi( t_n^ij )
            xt = fsw( t, se, eta )
!            
!           upper bound of switch region
            swthr = one + (se + 1.d0)*eta / 2.d0
!            
!           t_n^ij belongs to switch region
            if ( grad .and. ( t.lt.swthr .and. t.gt.swthr-eta ) ) then
!                    
              fac = dfsw( t, se, eta ) / rsph(jsph)
!    
!             accumulate for zi
!             -----------------
              zi(:,i,isph) = zi(:,i,isph) + fac*v(:)/vv
!              
            endif
!    
!           accumulate for fi
!           -----------------
            fi(i,isph) = fi(i,isph) + xt
!            
          enddo
!    
!         compute ui
!         ----------
          if ( fi(i,isph).le.one )  ui(i,isph) = one - fi(i,isph)
!
!
!---------------------------------------------------------------------------
! REMARK : the following correction would actually make zi(n,i) equal to
!          d ui(n,i) / dr_i , I think ...
!
!         if ( fi(i,isph).gt.one )  zi(:,i,isph) = zero
!---------------------------------------------------------------------------
!    
!    
        enddo
      enddo
!      
      !$omp end parallel do
!    
!    
!     build cavity array
!     ==================
!    
!     initialize number of cavity points
      ncav=0
!      
!     loop over spheres
      do isph = 1,nsph
!    
!       loop over integration points
        do i = 1,ngrid
!        
!         positive contribution from integration point
          if ( ui(i,isph).gt.zero ) then
!
!           accumulate
            ncav = ncav + 1
!                  
          endif
        enddo
      enddo
!    
!     allocate cavity array
      allocate( ccav(3,ncav) , stat=istatus )
      if ( istatus .ne. 0 ) then
        write(*,*)'ddinit : [3] allocation failed!'
        stop
      endif
!    
!     update memory usage
      memuse = memuse + 3*ncav
      memmax = max(memmax,memuse)
!
!     initialize cavity array index
      ii = 0
!
!     loop over spheres
      do isph = 1,nsph
!
!       loop over integration points
        do i = 1,ngrid
!
!         positive contribution from integration point
          if ( ui(i,isph).gt.zero ) then
!
!           advance cavity array index
            ii = ii + 1
!
!           store point
            ccav(:,ii) = csph(:,isph) + rsph(isph)*grid(:,i)
!            
          endif
        enddo
      enddo
!
!
!     safely initialize do_diag
!     =========================
      do_diag = .true.
!
!
!     derivatives d_j U_i^n
!     ---------------------
!
!     initialize
!
!          j n i
      du(:,:,:,:) = zero
!
!     loop over integration points
      do ig = 1, ngrid
!
!       loop over i-sphere
        do isph = 1,nsph
!   
!         loop over neighbors of i-sphere, i.e., potentially nonzero derivatives
          do ji = inl(isph), inl(isph+1) - 1
!          
!           neighbor is j-sphere, i.e., compute j-derivative
            jsph  = nl(ji)
!            
            vji   = csph(:,isph) + rsph(isph)*grid(:,ig) - csph(:,jsph)
            vvji  = sqrt(dot_product(vji,vji))
            tji   = vvji/rsph(jsph)
!            
!           switch region
            swthr = one + (se + 1.d0)*eta / 2.d0
!            
            if ( tji.lt.swthr .and. tji.gt.swthr-eta .and. fi(ig,isph).le.one ) then
!                    
              sji = vji/vvji
              fac = dfsw(tji,se,eta)/rsph(jsph)
!
!             accumulate
              du(1:3,jsph,ig,isph) = du(1:3,jsph,ig,isph) + fac*sji

              du(1:3,isph,ig,isph) = du(1:3,isph,ig,isph) - fac*sji

            endif
!            
          enddo
!
        enddo
      enddo

!!!      do ig = 1,ngrid
!!!        do isph = 1,nsph
!!!          do jsph = 1,nsph
!!!
!!!            write(*,3131) ig, isph, jsph, du(1:3,isph,ig,jsph) 
!!! 3131       format( ' n = ',i2,'; i = ',i2,'; j = ',i2,'; du = ',3(e12.5,2x) )
!!!
!!!          enddo
!!!        enddo
!!!      enddo
!!!      read(*,*)


!
      return
!
!
endsubroutine ddinit
!---------------------------------------------------------------------------------
!
!
!
!
!
!---------------------------------------------------------------------------------
subroutine memfree
!
      implicit none
      integer :: istatus, istatus0
!
!---------------------------------------------------------------------------------
!
!     initialize deallocation flags
      istatus0 = 0 ; istatus = 0
!      
!     deallocate the arrays
      if( allocated(rsph)   )  deallocate( rsph   , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(csph)   )  deallocate( csph   , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(ccav)   )  deallocate( ccav   , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(w)      )  deallocate( w      , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(grid)   )  deallocate( grid   , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(basis)  )  deallocate( basis  , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(inl)    )  deallocate( inl    , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(nl)     )  deallocate( nl     , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(fact)   )  deallocate( fact   , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(facl)   )  deallocate( facl   , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(facs)   )  deallocate( facs   , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(ui)     )  deallocate( ui     , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(du)     )  deallocate( du     , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(fi)     )  deallocate( fi     , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(zi)     )  deallocate( zi     , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(read_x) )  deallocate( read_x , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(read_y) )  deallocate( read_y , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(read_z) )  deallocate( read_z , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(read_r) )  deallocate( read_r , stat=istatus ) ; istatus0 = istatus0 + istatus
      if( allocated(read_q) )  deallocate( read_q , stat=istatus ) ; istatus0 = istatus0 + istatus
!
      if ( istatus0.ne.0 ) then
        write(*,*)'memfree : dallocation failed!'
        stop
      endif
!
!     update memory usage
      memuse = memuse - 4*nsph - 4*ngrid - nbasis*ngrid - nsph-1 - nsph*nngmax - &
               2*ngrid*nsph - 2*lmax-1 - 3*nbasis
      if ( grad )  memuse = memuse - 3*ngrid*nsph
!
!
endsubroutine memfree
!---------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------------------------
!
!
!
!
!
!------------------------------------------------------------------------------------------------
real*8 function fsw( t, s, eta )
!
      implicit none
      real*8, intent(in) :: t
      real*8, intent(in) :: s
      real*8, intent(in) :: eta
!
      real*8 :: a, b, flow, x
      real*8, parameter :: f6=6.0d0, f10=10.d0, f12=12.d0, f15=15.d0
!
!------------------------------------------------------------------------------------------------
!
!     shift :
!     s =  0   =>   t - eta/2  [ CENTERED ]
!     s =  1   =>   t - eta    [ EXTERIOR ]
!     s = -1   =>   t          [ INTERIOR ]
!
!     apply shift
      x = t - (s + 1.d0)*eta / 2.d0
!      
!     lower bound of switch region
      flow = one - eta
!      
!     define switch function \chi
      if     ( x.ge.one ) then
!              
        fsw = zero
!        
      elseif ( x.le.flow ) then
!      
        fsw = one
!        
      else
!              
        a = f15*eta - f12
        b = f10*eta*eta - f15*eta + f6
        fsw = ((x-one)*(x-one)*(one-x)*(f6*x*x + a*x + b))/(eta**5)
!        
      endif
!
!
endfunction fsw
!------------------------------------------------------------------------------------------------
!
!
!
!
!------------------------------------------------------------------------------------------------
real*8 function dfsw( t, s, eta )
!
      implicit none
      real*8, intent(in) :: t
      real*8, intent(in) :: s
      real*8, intent(in) :: eta
!      
      real*8  flow, x
      real*8, parameter :: f30=30.0d0
!
!------------------------------------------------------------------------------------------------
!
!     shift :
!     s =  0   =>   t - eta/2  [ CENTERED ]
!     s =  1   =>   t - eta    [ EXTERIOR ]
!     s = -1   =>   t          [ INTERIOR ]
!
!     apply shift
      x = t - (s + 1.d0)*eta / 2.d0
!
!     lower bound of switch region
      flow = one - eta
!      
!     define derivative of switch function \chi
      if     ( x.ge.one ) then
!              
        dfsw = zero
!              
      elseif ( x.le.flow ) then
!              
        dfsw = zero
!              
      else
!              
        dfsw = f30*(one-x)*(x-one)*(x-one+eta)*(x-one+eta)/(eta**5)
!              
      endif
!
!
endfunction dfsw
!------------------------------------------------------------------------------------------------
!
!
!
!
!
!
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
!
  subroutine ylmbas(x,basloc,vplm,vcos,vsin)
  implicit none
  real*8, dimension(3), intent(in) :: x
  real*8, dimension(nbasis), intent(out) :: basloc, vplm
  real*8, dimension(lmax+1), intent(out) :: vcos, vsin
  !
  integer :: l, m, ind
  real*8  :: cthe, sthe, cphi, sphi, plm

  basloc=zero ; vplm=zero ; vcos=zero ; vsin=zero
!
! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
! coordinates of x.
!
! evaluate cos( theta ) ; sin( theta )  
  cthe = x(3)
  sthe = sqrt(one - cthe*cthe)
!
! evalutate cos( phi ) ; sin( phi )
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
  if(sthe.ne.zero) then
    call trgev(cphi,sphi,vcos,vsin)
  else
    vcos = one
    vsin = zero
  end if
!
! evaluate the generalized legendre polynomials
  call polleg(cthe,sthe,vplm)
!
! now build the spherical harmonics. we will distinguish m=0,
! m>0 and m<0:
  do l = 0, lmax
    ind = l**2 + l + 1
!
!   m = 0
    basloc(ind) = facs(ind)*vplm(ind)
!
    do m = 1, l
!    
      plm = vplm(ind+m)
!
!     m > 0
      basloc(ind+m) = facs(ind+m)*plm*vcos(m+1)
!
!     m < 0
      basloc(ind-m) = facs(ind-m)*plm*vsin(m+1)
!
    end do
  end do
  return
  end subroutine ylmbas
  !
subroutine dbasis(x,basloc,dbsloc,vplm,vcos,vsin)
!
      implicit none
      real*8, dimension(3),        intent(in)    :: x
      real*8, dimension(nbasis),   intent(inout) :: basloc, vplm
      real*8, dimension(3,nbasis), intent(inout) :: dbsloc
      real*8, dimension(lmax+1),   intent(inout) :: vcos, vsin
!
      integer :: l, m, ind, VC, VS
      real*8  :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp
      real*8  :: et(3), ep(3)
!
!------------------------------------------------------------------------------------
!
!     get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
!     coordinates of x.
      cthe = x(3)
      sthe = sqrt(one - cthe*cthe)
!    
!     not ( NORTH or SOUTH pole )
      if ( sthe.ne.zero ) then
        cphi = x(1)/sthe
        sphi = x(2)/sthe
!    
!     NORTH or SOUTH pole
      else
        cphi = one
        sphi = zero
      end if
!
!     evaluate the derivatives of theta and phi:
!
      et(1) = cthe*cphi
      et(2) = cthe*sphi
      et(3) = -sthe

!     not ( NORTH or SOUTH pole )
      if( sthe.ne.zero ) then
        ep(1) = -sphi/sthe
        ep(2) = cphi/sthe
        ep(3) = zero
!        
!     NORTH or SOUTH pole
      else
        ep(1) = zero
        ep(2) = one
        ep(3) = zero
!        
      end if
!
!     evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
!     pointless if z = 1, as the only non vanishing terms will be the 
!     ones with m=0.
!
!     not ( NORTH or SOUTH pole )
      if ( sthe.ne.zero ) then
!              
        call trgev( cphi, sphi, vcos, vsin )
!        
!     NORTH or SOUTH pole
      else
!              
        vcos = one
        vsin = zero
!        
      end if
      VC=zero
      VS=cthe
!
!     evaluate the generalized legendre polynomials.
!
      call polleg( cthe, sthe, vplm )
!
!     now build the spherical harmonics. we will distinguish m=0,
!     m>0 and m<0:
!
      basloc = zero
      dbsloc = zero
      do l = 0, lmax
        ind = l*l + l + 1
      ! m = 0
        fln = facs(ind)   !  =>  fln (Fil) == vfac (Ben)
        basloc(ind) = fln*vplm(ind)
        if (l.gt.0) then
          dbsloc(:,ind) = fln*vplm(ind+1)*et(:)
        else
          dbsloc(:,ind) = zero
        end if
      !dir$ simd
        do m = 1, l
          fln = facs(ind+m)
          plm = fln*vplm(ind+m)   !  =>  plm (Fil) == fln (Fil) * plm (Ben) 
          pp1 = zero
          if (m.lt.l) pp1 = -pt5*vplm(ind+m+1)
          pm1 = pt5*(dble(l+m)*dble(l-m+1)*vplm(ind+m-1))
          pp  = pp1 + pm1   !  =>  pp (Fil) == dPlm (Ben)
!
!         m > 0
!         -----
!
          basloc(ind+m) = plm*vcos(m+1)
!          
!         not ( NORTH or SOUTH pole )
          if ( sthe.ne.zero ) then
! 
            dbsloc(:,ind+m) = -fln*pp*vcos(m+1)*et(:) - dble(m)*plm*vsin(m+1)*ep(:)
!
! Fil  =>                         fln           (-et              pp           vcos                         m ep               plm / fln   vsin                     )
! Ben  =>   grad( sind,:,ind+m) = vfac(1,ind+m)*(-ethe( sind,:).*(dPlm( sind).*vcos( sind,m+1)*ones(1,3)) - m*ephi( sind,:).*( plm( sind).*vsin(sind,m+1)*ones(1,3)));
!            
!         NORTH or SOUTH pole
          else
!                  
            dbsloc(:,ind+m) = -fln*pp*vcos(m+1)*et(:) - fln*pp*ep(:)*VC
!
! Fil =>                          fln           (-et              pp           vcos                           ep              pp                                    )
! Ben =>    grad(~sind,:,ind+m) = vfac(1,ind+m)*(-ethe(~sind,:).*(dPlm(~sind).*vcos(~sind,m+1)*ones(1,3)) -   ephi(~sind,:).*(dPlm(~sind).*            VC*ones(1,3)));
!
          endif
!
!         m < 0
!         -----
!
          basloc(ind-m) = plm*vsin(m+1)
!          
!         not ( NORTH or SOUTH pole )
          if ( sthe.ne.zero ) then
! 
            dbsloc(:,ind-m) = -fln*pp*vsin(m+1)*et(:) + dble(m)*plm*vcos(m+1)*ep(:)
!
!         NORTH or SOUTH pole
          else
!                  
            dbsloc(:,ind-m) = -fln*pp*vsin(m+1)*et(:) - fln*pp*ep(:)*VS
!            
          endif
!  
        end do
      end do
      return
!
!
endsubroutine dbasis
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
  real*8, dimension( max((lmax+1),2) ), intent(inout) :: cx, sx
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
!
! 
!---------------------------------------------------------------------
! Purpose : compute
!
!                       l
!     sum   4pi/(2l+1) t  * Y_l^m( s ) * sigma_l^m
!     l,m                                           
!
! which is need to compute action of COSMO matrix L.
!---------------------------------------------------------------------
!
real*8 function intmlp( t, sigma, basloc )
!  
      implicit none
      real*8, intent(in) :: t
      real*8, dimension(nbasis), intent(in) :: sigma, basloc
!
      integer :: l, ind
      real*8  :: tt, ss, fac
!
!---------------------------------------------------------------------
!
!     initialize t^l
      tt = one
!
!     initialize
      ss = zero
!
!     loop over l
      do l = 0,lmax
!      
        ind = l*l + l + 1
!
!       update factor 4pi / (2l+1) * t^l
        fac = tt / facl(ind)
!
!       contract over l,m and accumulate
        ss = ss + fac * dot_product( basloc(ind-l:ind+l), &
                                     sigma(ind-l:ind+l)   )
!
!       update t^l
        tt = tt*t
!        
      enddo
!      
!     redirect
      intmlp = ss
!
!
endfunction intmlp
!---------------------------------------------------------------------
!
!
!
!
subroutine wghpot( phi, g )
!
      implicit none
!
      real*8, dimension(ncav),       intent(in)  :: phi
      real*8, dimension(ngrid,nsph), intent(out) :: g
!
      integer isph, ig, ic
!
      ic = 0 ; g(:,:)=0.d0
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
!
!
endsubroutine wghpot
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
!
!
!-----------------------------------------------------------------------------------
! Purpose : compute
!
!   v_l^m = v_l^m +
!
!               4 pi           l
!     sum  sum  ---- ( t_n^ji )  Y_l^m( s_n^ji ) W_n^ji [ \xi_j ]_n
!      j    n   2l+1
!
! which is related to the action of the adjont COSMO matrix L^* in the following
! way. Set
!
!   [ \xi_j ]_n = sum  Y_l^m( s_n ) [ s_j ]_l^m
!                 l,m
!
! then
!
!   v_l^m = -   sum    (L^*)_ij s_j
!             j \ne i 
!
! The auxiliary quantity [ \xi_j ]_l^m needs to be computed explicitly.
!-----------------------------------------------------------------------------------
!
subroutine adjrhs1( isph, xi, vlm, basloc, vplm, vcos, vsin )
!
      implicit none
      integer,                       intent(in)    :: isph
      real*8, dimension(ngrid,nsph), intent(in)    :: xi
      real*8, dimension(nbasis),     intent(inout) :: vlm
      real*8, dimension(nbasis),     intent(inout) :: basloc, vplm
      real*8, dimension(lmax+1),     intent(inout) :: vcos, vsin
!
      integer :: ij, jsph, ig, l, ind, m
      real*8  :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t
!      
!-----------------------------------------------------------------------------------
!
!     loop over neighbors of i-sphere
      do ij = inl(isph),inl(isph+1)-1
!
!       j-sphere is neighbor
        jsph = nl(ij)
!
!       loop over integration points
        do ig = 1,ngrid
!        
!         compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
          vji  = csph(:,jsph) + rsph(jsph)*grid(:,ig) - csph(:,isph)
          vvji = sqrt(dot_product(vji,vji))
          tji  = vvji/rsph(isph)
!
!         point is INSIDE i-sphere (+ transition layer)
!         ---------------------------------------------
          if ( tji.lt.( one + (se+one)/two*eta ) ) then
!                  
!           compute s_n^ji
            sji = vji/vvji
!
!           compute \chi( t_n^ji )
            xji = fsw( tji, se, eta )
!
!           compute W_n^ji
            if ( fi(ig,jsph).gt.one ) then
!                    
              oji = xji/fi(ig,jsph)
!              
            else
!                    
              oji = xji
!              
            endif
!            
!           compute Y_l^m( s_n^ji )
            call ylmbas( sji, basloc, vplm, vcos, vsin )
!            
!           initialize ( t_n^ji )^l
            t = one
!            
!           compute w_n * xi(n,j) * W_n^ji
            fac = w(ig) * xi(ig,jsph) * oji
!            
!           loop over l
            do l = 0,lmax
!            
              ind  = l*l + l + 1
!
!             compute 4pi / (2l+1) * ( t_n^ji )^l * w_n * xi(n,j) * W_n^ji
              ffac = fac*t/facl(ind)
!
!             loop over m
              do m = -l,l
!              
                vlm(ind+m) = vlm(ind+m) + ffac*basloc(ind+m)
!                
              enddo
!
!             update ( t_n^ji )^l
              t = t*tji
!              
            enddo
!            
          endif
        enddo
      enddo
!
!
endsubroutine adjrhs1
!-----------------------------------------------------------------------------------

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
  1010 format( ' Parameters:',/, &
               '   number of grid points:                  '8x,i8,/,   &
               '   number of spheres:                      '8x,i8,/,   &
               '   lmax for the spherical harmonics basis: '8x,i8,/,   &
               '   convergence threshold:                  '8x,d8.1,/, &
               '   regularization parameters (eta,s):      ',f8.3,f8.3,/,&
               '   dielectric constant:                   ',e12.5/)
!               
  if ( iprint.gt.0 ) then
!          
    write(iout,1000)
    write(iout,1010) ngrid, nsph, lmax, 10.0d0**(-iconv), eta, se,eps
!    
    if ( iscrf.eq.0 ) then 
      write(iout,1011) 
 1011 format( ' Use COSMO. '/ )      
    else
      write(iout,1012) 
 1012 format( ' Use PCM. '/ )      
    endif
!
    if ( igrad.eq.1 )  write(iout,1013)
 1013 format( ' Compute forces.'// )   
!
  endif
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
  real*8  :: vvij, tij, xij, oij, t, fac, fl, f1, f2, f3, beta, tlow, thigh
  real*8  :: vij(3), sij(3), alp(3), va(3)
!
  tlow  = one - pt5*(one - se)*eta
  thigh = one + pt5*(one + se)*eta
!
  do ig = 1, ngrid
    va = zero
    do ij = inl(isph), inl(isph+1) - 1
      jsph = nl(ij)
      vij  = csph(:,isph) + rsph(isph)*grid(:,ig) - csph(:,jsph)
      vvij = sqrt(dot_product(vij,vij))
      tij  = vvij/rsph(jsph)
!
      if (tij.ge.thigh) cycle
!
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
      xij = fsw(tij,se,eta)
      if (fi(ig,isph).gt.one) then
        oij = xij/fi(ig,isph)
        f2  = -oij/fi(ig,isph)
      else
        oij = xij
        f2  = zero
      end if
      f1 = oij/rsph(jsph)
      va(:) = va(:) + f1*alp(:) + beta*f2*zi(:,ig,isph)
      if (tij .gt. tlow) then
!!!      if (tij .gt. (one-eta*rsph(jsph))) then
!!!        f3 = beta*dfsw(tij,eta*rsph(jsph))/rsph(jsph)
        f3 = beta*dfsw(tij,se,eta)/rsph(jsph)
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
  real*8  :: vvji, tji, xji, oji, t, fac, fl, f1, f2, beta, di, tlow, thigh
  real*8  :: b, g1, g2, vvjk, tjk, f, xjk
  real*8  :: vji(3), sji(3), alp(3), vb(3), vjk(3), sjk(3), vc(3)
!
  tlow  = one - pt5*(one - se)*eta
  thigh = one + pt5*(one + se)*eta
!
  do ig = 1, ngrid
    vb = zero
    vc = zero
    do ji = inl(isph), inl(isph+1) - 1
      jsph = nl(ji)
      vji  = csph(:,jsph) + rsph(jsph)*grid(:,ig) - csph(:,isph)
      vvji = sqrt(dot_product(vji,vji))
      tji  = vvji/rsph(isph)
!
      if (tji.gt.thigh) cycle
!
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
      xji = fsw(tji,se,eta)
      if (fi(ig,jsph).gt.one) then
        oji = xji/fi(ig,jsph)
      else
        oji = xji
      end if
      f1 = oji/rsph(isph)
      vb = vb + f1*alp*xi(ig,jsph)
      if (tji .gt. tlow) then
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
              if (tjk .le. thigh) then
                proc = .true.
                sjk  = vjk/vvjk
                call ylmbas(sjk,basloc,vplm,vcos,vsin)
                g1  = intmlp(tjk,sigma(:,ksph),basloc)
                xjk = fsw(tjk,se,eta)
                b   = b + g1*xjk
              end if
            end if
          end do
          if (proc) then
!!!            g1 = di*di*dfsw(tji,eta*rsph(isph))/rsph(isph)
            g1 = di*di*dfsw(tji,se,eta)/rsph(isph)
            g2 = g1*xi(ig,jsph)*b
            vc = vc + g2*sji
          end if
        else
          di  = one
          fac = zero
        end if
!!!        f2 = (one-fac)*di*dfsw(tji,eta*rsph(isph))/rsph(isph)
        f2 = (one-fac)*di*dfsw(tji,se,eta)/rsph(isph)
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
      swthr = one + (se + 1.d0)*eta / 2.d0
      if (tji.lt.swthr .and. tji.gt.swthr-eta .and. ui(ig,jsph).gt.zero) then
        sji = vji/vvji
!!!        fac = - dfsw(tji,eta*rsph(isph))/rsph(isph)
        fac = - dfsw(tji,se,eta)/rsph(isph)
        alp = alp + fac*phi(ig,jsph)*xi(ig,jsph)*sji
      end if
    end do
    fx = fx - w(ig)*alp
  end do
  return 
  end subroutine fdoga
!
!
!-------------------------------------------------------------------------------
!
!
!
!
!---------------------------------------------------------------------
! Purpose : compute
!                                       
!     sum   4pi/(2l'+1) * Y_l'^m'(s_n) * sigma_l'^m'
!    l',m'                             
!
!---------------------------------------------------------------------
real*8 function cstmlp( sigma, basloc )
!
      implicit none
      real*8, dimension(nbasis), intent(in) :: sigma, basloc
!
      real*8  :: ss
!
!---------------------------------------------------------------------
!
      ss = dot_product( basloc(:)/facl(:), sigma(:) )
      cstmlp = ss
!      
      return
!      
endfunction cstmlp
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
! Purpose : compute :
!
!        4 pi (l+1)
!   sum  ----------  t^l  \nu_l^m  Y_l^m
!   l,m    2l + 1
!
! Remark : this expression is used for EXTENSION of potential
!          outside sphere. Index l starts from 0 !
!---------------------------------------------------------------------
real*8 function dtslm( t, nu, basloc )
!
      implicit none
      real*8,                    intent(in) :: t
      real*8, dimension(nbasis), intent(in) :: nu
      real*8, dimension(nbasis), intent(in) :: basloc
!
      integer :: l, ind, m
      real*8  :: fl, fac, ss, tt
!
!---------------------------------------------------------------------
!
!     initialize
      ss = zero ; tt = one
!
!     loop over degree of spherical harmonics 
      do l = 0, lmax
!      
!       index associated to Y_l^0
        ind = l*l + l + 1
!        
!       compute factor
        fl  = dble(l)
        fac = four*pi*(fl+one)/(two*fl + one)*tt
!
!       loop over order of spherical harmonics
        do m = -l, l
!
!         accumulate
          ss = ss + fac*nu(ind+m)*basloc(ind+m)
!          
        end do
!
!       compute t^l
        tt = tt*t
!        
      end do
!      
!     return value
      dtslm = ss
!
!
endfunction dtslm
!---------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------
! Purpose : compute
!
!        4\pi l
!   sum  ------ 1/t^(l+1) \nu_l^m Y_l^m
!   l,m  2l + 1
!
! which is needed to compute the action of PCM .
!---------------------------------------------------------------------
!
real*8 function dtslm2( t, nu, basloc )
!
      implicit none
      real*8,                    intent(in) :: t                   
      real*8, dimension(nbasis), intent(in) :: nu
      real*8, dimension(nbasis), intent(in) :: basloc
!
      integer :: l, ind, m
      real*8  :: fl, fac, ss, tt
!      
!---------------------------------------------------------------------
!
!     initialize 1/t^(l+1)
      tt = one/t
!      
!     initialize return value
      ss = zero
!
!     loop over degree of spherical harmonics 
      do l = 0,lmax
!
        ind = l*l + l + 1
!
!       compute factor
        fl  = dble(l)
        fac = four*pi*fl/(two*fl + one)*tt
!
!       loop over order of spherical harmonics        
        do m = -l,l
!
!         accumulate over l, m
          ss = ss + fac * nu(ind+m) * basloc(ind+m)
!          
        enddo
!
!       update 1/t^(l+1)
        tt = tt/t
!        
      enddo
!
!     return value
      dtslm2 = ss
!
!
endfunction dtslm2
!---------------------------------------------------------------------
!
!
!
!
!---------------------------------------------------------------------
! Purpose : compute
!
!          4\pi l
!   xlm =  ------ 1/t^(l+1) Y_l^m
!          2l + 1
!
! which is needed to compute the action of adjoint PCM .
!---------------------------------------------------------------------
!
subroutine tylm( t, basloc, xlm )
!
      implicit none
      real*8,                    intent(in ) :: t                   
      real*8, dimension(nbasis), intent(in ) :: basloc
      real*8, dimension(nbasis), intent(out) :: xlm
!
      integer :: l, ind, m
      real*8  :: fl, fac, ss, tt
!      
!---------------------------------------------------------------------
!
!     initialize 1/t^(l+1)
      tt = one/t
!      
!     initialize
      xlm = zero
!
!     loop over degree of spherical harmonics 
      do l = 0,lmax
!
        ind = l*l + l + 1
!
!       compute factor
        fl  = dble(l)
        fac = four*pi*fl/(two*fl + one)*tt
!
!       loop over order of spherical harmonics        
        do m = -l,l
!
!         accumulate over l, m
          xlm(ind+m) = fac * basloc(ind+m)
!          
        enddo
!
!       update 1/t^(l+1)
        tt = tt/t
!        
      enddo
!
!
endsubroutine tylm
!---------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
! Purpose : compute
!
!   \Phi( n ) =
!     
!                       4 pi           l
!     sum  W_n^ij  sum  ---- ( t_n^ij )  Y_l^m( s_n^ij ) [ \sigma_j ]_l^m
!      j           l,m  2l+1
!
! which is related to the action of the COSMO matrix L in the following
! way :
!
!   -   sum    L_ij \sigma_j = sum  w_n Y_l^m( s_n ) \Phi( n ) 
!     j \ne i                   n
!
! This second step is performed by routine "intrhs".
!------------------------------------------------------------------------
!
subroutine calcv2( first, isph, pot, sigma, basloc, vplm, vcos, vsin )
!
      logical,                        intent(in)    :: first
      integer,                        intent(in)    :: isph
      real*8, dimension(nbasis,nsph), intent(in)    :: sigma
      real*8, dimension(ngrid),       intent(inout) :: pot
      real*8, dimension(nbasis),      intent(inout) :: basloc
      real*8, dimension(nbasis),      intent(inout) :: vplm
      real*8, dimension(lmax+1),      intent(inout) :: vcos
      real*8, dimension(lmax+1),      intent(inout) :: vsin
!
      integer :: its, ij, jsph
      real*8  :: vij(3), sij(3)
      real*8  :: vvij, tij, xij, oij, stslm, stslm2, stslm3
!
!------------------------------------------------------------------------
!
!     initialize
      pot(:) = zero
!
!     if 1st iteration of Jacobi method, then done!
      if ( first )  return
!
!     loop over grid points
      do its = 1,ngrid
!
!       contribution from integration point present
        if ( ui(its,isph).lt.one ) then
!
!         loop over neighbors of i-sphere
          do ij = inl(isph),inl(isph+1)-1
!
!           neighbor is j-sphere
            jsph = nl(ij)
!            
!           compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
            vij  = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
            vvij = sqrt( dot_product( vij, vij ) )
            tij  = vvij / rsph(jsph) 
!
!           compute s_n^ij = ( r_i + \rho_i s_n - r_j ) / | ... |
            sij = vij / vvij
!            
!           compute \chi( t_n^ij )
            xij = fsw( tij, se, eta )
!
!           compute W_n^ij
            if ( fi(its,isph).gt.one ) then
!
              oij = xij / fi(its,isph)
!
            else
!
              oij = xij
!
            endif
!
!           point is INSIDE j-sphere
!           ------------------------
            if ( tij.lt.one ) then
!
!             compute Y_l^m( s_n^ij )
              call ylmbas( sij, basloc, vplm, vcos, vsin )
!                    
!             accumulate over j, l, m
              pot(its) = pot(its) + oij * intmlp( tij, sigma(:,jsph), basloc )
!              
!           point is OUTSIDE j-sphere (+ transition layer) [EXTENSION]
!           ----------------------------------------------------------
            elseif ( tij.lt.( one + (se+one)/two*eta ) ) then
!                    
!             extension of potential
              select case(ext0)
!
!             t^l extension
              case(0)
!
!             compute Y_l^m( s_n^ij )
              call ylmbas( sij, basloc, vplm, vcos, vsin )
!              
!             accumulate over j, l, m
              pot(its) = pot(its) + oij*intmlp(tij,sigma(:,jsph),basloc)
!
!             constant extension
              case(1)
              pot(its) = pot(its) + oij*cstmlp(    sigma(:,jsph),basloc)
!
              endselect
!
!                    
!!!            elseif ( tij.lt.( one + (se+one)/two*eta ) ) then
!!!!                    
!!!!             compute l,m-th row of L_jj * X_j vector
!!!              pot(its) = pot(its) + oij * cstmlp(      sigma(:,jsph), basloc )
!!!!                    
!!!      !       pot(its) = pot(its) + oij*extmlp(tij,sigma(:,jsph),basloc)
!!!      !       pot(its) = pot(its) + oij*stslm2(lmax,tij,sigma(:,jsph),basloc)
!!!      !       pot(its) = pot(its) + oij*stslm3(lmax,tij,sigma(:,jsph),basloc)
!!!!                    
            endif
          end do
        end if
      end do
!      
      return
!      
!      
endsubroutine calcv2
!
!
endmodule ddcosmo
