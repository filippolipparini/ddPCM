subroutine debug_tests()
!
      use ddcosmo , only : lmax, ngrid, nsph, read_x, read_y, read_z, read_r, &
                           ddinit, ncav, nbasis, read_q, ccav, memfree, &
                           reset_ngrid, iquiet, iprint, args, &
                           read_molecule_file
!      
      implicit none
      real*8 :: xx(1), esolv, charge(nsph)!, x(nsph), y(nsph), z(nsph), rvdw(nsph)
      real*8, allocatable :: phi(:), psi(:,:), f(:,:), sigma(:,:), s(:,:), &
                             phi_eps(:,:),grid_aux(:,:),w_aux(:)
      integer :: idec, igrid, lmax_save, ngrid_save, istatus, i, nsph_save, &
                 iprint_save, np, iidec
      logical :: iquiet_save
!      
      integer, parameter, dimension(32) :: ngrid_vec = (/   6,  14,  26,  38,  50,  74,  86, 110,  &
                                                          146, 170, 194, 230, 266, 302, 350, 434,  &
                                                          590, 770, 974,1202,1454,1730,2030,2354,  &
                                                         2702,3074,3470,3890,4334,4802,5294,5810/)
!
!----------------------------------------------------------------------------------
!
!     set on quiet mode
      iquiet_save = iquiet ; iquiet = .true.
      iprint_save = iprint ; iprint = 0
!      
!     display option menu in infinite loop
      do
!
        write(*,*)'==========================================================='
        write(*,*)'TESTS MENU :                                               '
        write(*,*)''
        write(*,*)'Check derivatives of switch function U_j .................1'
        write(*,*)'Check derivatives of COSMO matrix L ......................2'
        write(*,*)'Check derivatives of PCM   matrix A ......................3'
        write(*,*)'Check COSMO forces over angular momenta and grids ........4'
        write(*,*)'Check PCM   forces over angular momenta and grids ........5'
        write(*,*)'Count points in switch region ............................6'
        write(*,*)'Action of PCM with FMM ...................................7'
        write(*,*)'Compare PCM implementations ..............................8'
        write(*,*)'Check adjoint PCM ........................................9'
        write(*,*)''
        write(*,*)'Exit .....................................................0'
        write(*,*)'==========================================================='
!        
        read(*,*) idec
!
        select case(idec)
!
        case(0) ; exit
        case(1) ; call check_der_U()
        case(2) ; call check_der_L()
        case(3) ; call check_der_A()
!                
!       COSMO forces
!       ============
        case(4)
!
!         save parameters
          lmax_save = lmax ; ngrid_save = ngrid 
!
!         loop over angular momenta
          do lmax = 2,10
!
!           set ngrid so that 2*lmax is integrated exactely
            call reset_ngrid( igrid )
!
!           loop over extra grids
            do i = 1,4
!
!             free the memory
              call memfree()
!
!             recreate nsph, read_x, read_y, read_z, read_r, read_q
              call read_molecule_file()
!
!             redirect charge [ THIS WILL NEED TO BE CLEANED UP ]
              charge = read_q
!
!             initialize datastructure
              call ddinit( nsph, read_x, read_y, read_z, read_r )
!
!             allocate workspaces
              allocate( phi(ncav), psi(nbasis,nsph), sigma(nbasis,nsph), s(nbasis,nsph), f(3,nsph) , stat=istatus )
              if ( istatus.ne.0 ) then
                write(*,*)'debug_tests : COSMO allocation failed !'
                stop
              endif
!              
!             solute electrostatic potential phi, and psi vector
              call mkrhs( nsph, charge, read_x, read_y, read_z, ncav, ccav, phi, nbasis, psi )
!
!             cosmo equation
              call cosmo( .false., .true., phi, xx, psi, sigma, esolv )
!     
!             cosmo adjoint equation
              call cosmo( .true., .false., phi, xx, psi, s, esolv )
!              
!             compute forces
              call forces( nsph, charge, phi, sigma, s, f )
!
!             FD convergence test
              call check_forcesCOSMO( esolv, charge, f )
!              
!             deallocate workspaces
              deallocate( phi, psi, sigma, s, f , stat=istatus )
              if ( istatus.ne.0 ) then
                write(*,*)'debug_tests : COSMO deallocation failed !'
                stop
              endif
!
!             increment grid number        
              igrid = igrid + 1
! 
!             update number of grid points
              ngrid = ngrid_vec(igrid)
!              
            enddo
          enddo
!          
!         free the memory
          call memfree()
!
!         restore parameters
          lmax = lmax_save ; ngrid = ngrid_save
!          
!         recreate nsph, read_x, read_y, read_z, read_r, read_q
          call read_molecule_file()
!          
!         initialize datastructure
          call ddinit( nsph, read_x, read_y, read_z, read_r )
!                
!       PCM forces
!       ==========
        case(5)
!                
!         save parameters
          lmax_save = lmax ; ngrid_save = ngrid
!          
!         loop over angular momenta
          do lmax = 2,10
!
!           set ngrid so that 2*lmax is integrated exactely
            call reset_ngrid( igrid )
!
!           loop over extra grids
            do i = 1,4
!
!             free the memory
              call memfree()
!              
!             recreate nsph, read_x, read_y, read_z, read_r, read_q
              call read_molecule_file()
!
!             redirect charge [ THIS WILL NEED TO BE CLEANED UP ]
              charge = read_q
!
!             initialize datastructure
              call ddinit( nsph, read_x, read_y, read_z, read_r )
!
!             allocate workspaces
              allocate( phi(ncav), psi(nbasis,nsph), sigma(nbasis,nsph), s(nbasis,nsph), &
                        f(3,nsph) , phi_eps(nbasis,nsph) , stat=istatus )
              if ( istatus.ne.0 ) then
                write(*,*)'debug_tests : PCM allocation failed !'
                stop
              endif
!              
!             solute electrostatic potential phi, and psi vector
              call mkrhs( nsph, charge, read_x, read_y, read_z, ncav, ccav, phi, nbasis, psi )
!
!             pcm equations
              call pcm(   .false.,  .true., .true., phi,      xx, phi_eps )
              call cosmo( .false., .false.,          xx, phi_eps, psi, sigma, esolv )
!     
!             cosmo adjoint equation
              call cosmo( .true., .false., phi, xx, psi, s, esolv )
!              
!             compute forces
              call compute_forces( phi, charge, psi, sigma, phi_eps, f )
!
!             FD convergence test
              call check_forcesPCM( charge, f, esolv )
!              
!             deallocate workspaces
              deallocate( phi, psi, sigma, s, f, phi_eps , stat=istatus )
              if ( istatus.ne.0 ) then
                write(*,*)'debug_tests : PCM deallocation failed !'
                stop
              endif
!
!             increment grid number        
              igrid = igrid + 1
! 
!             update number of grid points
              ngrid = ngrid_vec(igrid)
!              
            enddo
          enddo
!          
!         free the memory
          call memfree()
!
!         restore parameters
          lmax = lmax_save ; ngrid = ngrid_save 
!          
!         recreate nsph, read_x, read_y, read_z, read_r, read_q
          call read_molecule_file()
!          
!         initialize datastructure
          call ddinit( nsph, read_x, read_y, read_z, read_r )
!
!       POINTS IN SWITCH REGION
!       =======================
        case(6)
!
!         header
          write(*,*)''
          write(*,*)'-----------------------'
          write(*,*)' lmax | ngrid | points '
          write(*,*)'-----------------------'
!                
!         save parameters
          lmax_save = lmax ; ngrid_save = ngrid
!          
!         loop over angular momenta
          do lmax = 2,10
!
!           set ngrid so that 2*lmax is integrated exactely
            call reset_ngrid( igrid )
!
!           loop over extra grids
            do i = 1,4
!
!             free the memory
              call memfree()
!              
!             recreate nsph, read_x, read_y, read_z, read_r, read_q
              call read_molecule_file()
!
!             redirect charge [ THIS WILL NEED TO BE CLEANED UP ]
              charge = read_q
!
!             initialize datastructure
              call ddinit( nsph, read_x, read_y, read_z, read_r )
!
!             count points
              call npoints_switch_region( np )
              write(*,1000) lmax, ngrid, np
 1000         format(2x,i4,3x,i5,3x,i6)         
!
!             increment grid number        
              igrid = igrid + 1
! 
!             update number of grid points
              ngrid = ngrid_vec(igrid)
!              
            enddo
          enddo
!
!         footer
          write(*,*)'-----------------------'
          write(*,*)''
!
!         repeat in infinite loop
          do 
!
            write(*,*)'Print points in Lebedev grid? ngrid - Yes ; 0 - No'
            read(*,*) iidec
!
!           exit inifinite loop
            if ( iidec.le.0 )  exit
!
            allocate( w_aux(iidec), grid_aux(3,iidec) , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'debug_tests : GRID failed allocation !'      
              stop
            endif
!
!           load grid
            call llgrid( iidec, w_aux, grid_aux )
!
!           print grid
            do i = 1,iidec
!
              write(*,1001) i,grid_aux(:,i)
 1001         format(' i = ',i3,' ; point = ',3(f6.3,2x) )            
!
            enddo
            write(*,*)''
!            
            deallocate( w_aux, grid_aux , stat=istatus )
            if ( istatus.ne.0 ) then
              write(*,*)'debug_tests : GRID failed deallocation !'      
              stop
            endif
!
          enddo
!          
!         free the memory
          call memfree()
!
!         restore parameters
          lmax = lmax_save ; ngrid = ngrid_save 
!          
!         recreate nsph, read_x, read_y, read_z, read_r, read_q
          call read_molecule_file()
!          
!         initialize datastructure
          call ddinit( nsph, read_x, read_y, read_z, read_r )
!
!       FMM
!       ===
        case(7) ; call test_fmm()
!
!       OLD AND NEW PCM
!       ===============
        case(8) ; call check_old_new_pcm()
!
!       CHECK ADJOINT PCM MATRIX
!       ========================
        case(9) ; call ADJcheck()
!          
        endselect
!                
      enddo
!
!     restore iquiet and iprint
      iquiet = iquiet_save
      iprint = iprint_save
!      
!
endsubroutine debug_tests
!----------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------
! Purpose : check that adjoint PCM matrix is formed correctly , and solutin of
!           adjoint problem is consistent with solution of direct problem.
!----------------------------------------------------------------------------------
!
subroutine ADJcheck()
!       
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, zero, iquiet, eps, iprint, &
                           one, do_diag
!
      implicit none
      real*8 :: f(nbasis,nsph), Af(nbasis), e(nbasis,nsph), s(nbasis,nsph)
      real*8 :: A(nbasis*nsph,nbasis*nsph), AT(nbasis*nsph,nbasis*nsph)
      real*8 :: xlm(nbasis), x(ngrid), vplm(nbasis), vcos(lmax+1), vsin(lmax+1), &
                basloc(nbasis), rvoid(1), s1, s2
      integer :: isph,jsph, i, j, ibeg, iend, iprint_save
      logical :: iquiet_save, do_diag_save
!
!--------------------------------------------------------------------------------
!
!     store flags
      iquiet_save = iquiet ; iprint_save = iprint ; do_diag_save = do_diag
!
!     activate flags
      iquiet = .true. ; iprint = 0 ; do_diag = .true.
!
!
!     check A = (A^T)^T
!     -----------------
!
!     construct A, A^T
      do isph = 1,nsph
!
        do jsph = 1,nsph
          do j = 1,nbasis
!
!           standard basis vector e_j
            e(:,   :) = zero
            e(j,jsph) = one
!
!           compute A  _i e_j
            ibeg = (isph-1)*nbasis+1
            iend = (isph-1)*nbasis+nbasis
            call mkrvec( isph, eps, e, A( ibeg:iend,(jsph-1)*nbasis+j), xlm, x, basloc, vplm, vcos, vsin )
!
!           compute A^T_i e_j
            call ADJvec( isph, eps, e, AT(ibeg:iend,(jsph-1)*nbasis+j), xlm, x, basloc, vplm, vcos, vsin )
!        
          enddo
        enddo
      enddo
!
!     initialize
      s1 = zero ; s2 = zero
!      
!     compute Frobenious norm of A and A^T
      do isph = 1,nsph
        do i = 1,nbasis
          do jsph = 1,nsph
            do j = 1,nbasis
!
!             accumulate A_ij^2
              s1 = s1 +   A ( (isph-1)*nbasis+i , (jsph-1)*nbasis+j )**2

!             accumulate ( A_ij - (A^T)_ij )^2
              s2 = s2 + ( A ( (isph-1)*nbasis+i , (jsph-1)*nbasis+j ) -  &
                          AT( (jsph-1)*nbasis+j , (isph-1)*nbasis+i ) )**2
!
            enddo
          enddo
        enddo
      enddo
!
!     the much neglected square roots 
      s1 = sqrt(s1) ; s2 = sqrt(s2)
!      
!     print
      write(*,1002) s2 / s1
 1002 format(' ||A - (A^T)^T||_F | / ||A||_F = ',e12.5)
      write(*,*) ''
!
!
!     extra print out
      if ( (s2/s1).ge.1.E-12 ) then
!   
!       clean up machine zeros
        do isph = 1,nsph
          do jsph = 1,nsph
            do i = 1,nbasis
              do j = 1,nbasis
!
                if ( abs(A(  (isph-1)*nbasis+i , (jsph-1)*nbasis+j )) .lt. 1.E-12 ) then
                  A(  (isph-1)*nbasis+i , (jsph-1)*nbasis+j ) = zero
                endif

                if ( abs(AT( (isph-1)*nbasis+i , (jsph-1)*nbasis+j )) .lt. 1.E-12 ) then
                  AT( (isph-1)*nbasis+i , (jsph-1)*nbasis+j ) = zero
                endif
               
              enddo 
            enddo 
          enddo 
        enddo 
!
!       print A
        write(*,*) 'A = '
        do isph = 1,nsph
          do i = 1,nbasis
            write(*,"(4x,300(e12.5,2x))") ( A((isph-1)*nbasis+i,j), j=1,nbasis*nsph )
          enddo
        enddo
        write(*,*)''
!
!       print A^T
        write(*,*) '(A^T)^T = '
        do isph = 1,nsph
          do i = 1,nbasis
            write(*,"(4x,300(e12.5,2x))") ( AT(j,(isph-1)*nbasis+i), j=1,nbasis*nsph )
          enddo
        enddo
        write(*,*)''
!        
      endif
!
!
!     check solution of adjoint problem : < e , f > = < A^T s , f > = < s , A f >
!     ---------------------------------------------------------------------------
!
!     initialize random number generator
      call random_seed
!      
!     build e, f 
      s1=zero ; s2=zero
      do isph = 1,nsph
        do j = 1,nbasis
!        
          call random_number( e(j,isph) )
          call random_number( f(j,isph) )
!          
          s1 = s1 + e(j,isph)**2
          s2 = s2 + f(j,isph)**2
!          
        enddo
      enddo
!
!     normalize
      e(:,:)=e(:,:)/sqrt(s1)
      f(:,:)=f(:,:)/sqrt(s2)
!
!     solve A_eps^T s = e
      call pcm( .true., .false., .true., rvoid, e, s )
!      
!     initialize
      Af(:) = zero ; s1 = zero ; s2 = zero
!      
      do isph = 1,nsph
!
!       compute A_i f 
        call mkrvec( isph, eps, f, Af(:), xlm, x, basloc, vplm, vcos, vsin )
!
!       accumulate < s , A f >
        s1 = s1 + dot_product( s(:,isph), Af(:) )
!        
!       accumulate < e , f >
        s2 = s2 + dot_product( e(:,isph), f(:,isph) )
!        
      enddo
!
!     print
      write(*,1001) abs(s2-s1) / abs(s1)
 1001 format(' | <e,f> - <s,Af> | / |<s,Af>| = ', e12.5)
      write(*,*) ''
!
!     restore flags
      iquiet = iquiet_save ; iprint = iprint_save ; do_diag = do_diag_save
!
!
endsubroutine ADJcheck
!---------------------------------------------------------------------------------------
!
!
!
!
!
!
!---------------------------------------------------------------------------------------
! Purpose : check rates of convergence of COSMO forces.
!---------------------------------------------------------------------------------------
subroutine check_forcesCOSMO( E0, charge, f )
!
      use ddcosmo , only : nbasis, nsph, iquiet, csph, rsph, memfree, ddinit, &
                           ncav, ccav, ngrid, zero, one, lmax, &
                           iprint, tokcal
!                           
      implicit none
      real*8,                    intent(in) :: E0
      real*8, dimension(  nsph), intent(in) :: charge
      real*8, dimension(3,nsph), intent(in) :: f
!
      integer,parameter :: niter = 4
!
      real*8 :: phi(ncav), psi(nbasis,nsph), sigma(nbasis,nsph), xx(1)
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), phi_eps(nbasis,nsph)
      real*8 :: x(nsph), y(nsph), z(nsph), rwork(niter,nsph*3),rrate(niter,nsph*3) , &
                hwork(niter,nsph*3),ework(niter,nsph*3)
      real*8 :: E_plus, err, eeps, h
      integer :: iter, icomp, ksph, nsph_save, j, ncav_save, nbasis_save
!      
      character(len=10) :: x1,x2
      character(len=30) :: fname
      integer :: fp, iprint_save, np
      logical :: iquiet_save
!
!---------------------------------------------------------------------------------------
!
!     store flags
      iquiet_save = iquiet ; iprint_save = iprint
!      
!     activate quiet flag
      iquiet = .true. ; iprint = 0
!
!     compute number of points in switch region
      call npoints_switch_region( np )
!
!     save initial DS
      nsph_save = nsph
      ncav_save = ncav
      nbasis_save = nbasis
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
!     set initial increment
      eeps=0.001d0
!
!     initialize
      rwork = zero ; rrate = zero ; ework = zero ; hwork = zero
!
!     loop over increments
      do iter = 1,niter
!        
!       loop over d / dr_k
        do ksph = 1,nsph_save
!
!         loop over components of d / dr_k
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph)*(1.d0+eeps) ; h = eeps*x_save(ksph)
            case(2) ; y(ksph) = y_save(ksph)*(1.d0+eeps) ; h = eeps*y_save(ksph)
            case(3) ; z(ksph) = z_save(ksph)*(1.d0+eeps) ; h = eeps*z_save(ksph)
            endselect
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!            
!           potential Phi and Psi vector
            call mkrhs( nsph_save, charge, x, y, z, ncav_save, ccav, phi, nbasis_save, psi )
!
!           solve COSMO equations       
            E_plus = zero ; sigma = zero
            call cosmo( .false., .true., phi, xx, psi, sigma, E_plus )
!
!           account for a potentially null shift
            if ( abs( h ).gt.1.E-12 ) then
!                    
!             numerical derivative of energy
              ework(iter,(ksph-1)*3+icomp) = ( E_plus - E0 ) / h
!              
!             account for a potentially null component of the force
              if ( abs( f(icomp,ksph) ).gt.1.E-12 ) then
!                      
!               compute relative error
                err = abs( (E_plus - E0) / h + f(icomp,ksph) ) / abs( f(icomp,ksph) )
!
!               store
                rwork(iter,(ksph-1)*3+icomp) = err
                hwork(iter,(ksph-1)*3+icomp) = h
!                
!               compute rate
                if ( iter.gt.1 ) then 
                  rrate(iter,(ksph-1)*3+icomp) =  log( rwork(iter-1,(ksph-1)*3+icomp) / &
                                                       rwork(iter  ,(ksph-1)*3+icomp)   ) / &
                                                  log( hwork(iter  ,(ksph-1)*3+icomp) / & 
                                                       hwork(iter-1,(ksph-1)*3+icomp)   )
                endif
              endif
            endif
!
          enddo
        enddo
!
        eeps = eeps / 2.d0
!
      enddo
!      
!     print header
      write(*,*)'-----------------------------------------------------------------------------------------------'
      if ( np.gt.0 ) then 
        write(*,7324) lmax,ngrid
      else
        write(*,7325) lmax,ngrid
      endif
 7324 format(' lmax = ',i2,' , ngrid = ',i4,' - Integration points PRESENT in the switch region')     
 7325 format(' lmax = ',i2,' , ngrid = ',i4,' - Integration points NOT PRESENT in the switch region')     
      write(*,*)''
!
!     print numerical derivative of energy
      write(*,*)'Numerical derivative of energy : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( ework(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!     print relative error
      write(*,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!     print rate of convergence
      write(*,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''

      write(*,1000)
 1000 format(1x,'ddCOSMO forces (atomic units):',/, &
                1x,' atom',15x,'x',15x,'y',15x,'z')
!                
      do ksph = 1, nsph
        write(*,'(1x,i5,3f16.8)') ksph, f(:,ksph)
      enddo
!
      write (*,'(1x,a,f14.6)') 'ddcosmo electrostatic solvation energy (kcal/mol):', E0*tokcal
      write(*,*) ''
!
!     print to file
!     -------------
!
      write(x1,'(I2.2)') lmax
      write(x2,'(I4.4)') ngrid
      fname = 'lmax' // trim(x1) // '_ngrid' // trim(x2)
      fp    = 17
!
      open( unit=fp, file=fname, form='formatted', access='sequential', status='unknown')

!     print header
      if ( np.gt.0 ) then 
        write(fp,7324) lmax,ngrid
      else
        write(fp,7325) lmax,ngrid
      endif
      write(fp,*)''
!
!      
      write(fp,*)'lmax,ngrid = ',lmax,ngrid
      write(fp,*) ''
!
!     print numerical derivative of energy
      write(fp,*)'Numerical derivative of energy : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( ework(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!
!     print relative error
      write(fp,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!      
!     print rate of convergence
      write(fp,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
!
      close(fp)
!      
!     restore DS
      call memfree
      call ddinit( nsph_save, x_save, y_save, z_save, r_save )
!
!     restore flags
      iquiet = iquiet_save ; iprint = iprint_save
!
!
endsubroutine check_forcesCOSMO
!---------------------------------------------------------------------------------------
!
!
!
!
!
!---------------------------------------------------------------------------------------
! Purpose : check rates of convergence of PCM forces.
!---------------------------------------------------------------------------------------
!
subroutine check_forcesPCM( charge, f, esolv )
!
      use ddcosmo , only : nbasis, nsph, iquiet, csph, rsph, memfree, ddinit, &
                           eps, ncav, ccav, ngrid, zero, sprod, wghpot, one, &
                           lmax, iprint
!                           
      implicit none
      real*8, dimension(  nsph), intent(in) :: charge
      real*8, dimension(3,nsph), intent(in) :: f
      real*8                                :: esolv
!
      integer,parameter :: niter = 4
!
      real*8 :: phi(ncav), psi(nbasis,nsph), g(ngrid,nsph), sigma(nbasis,nsph), &
                xx(1)
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), phi_eps(nbasis,nsph)
      real*8 :: x(nsph), y(nsph), z(nsph), rwork(niter,nsph*3), rrate(niter,nsph*3), &
                hwork(niter,nsph*3),ework(niter,nsph*3),fwork(niter,nsph*3),eework(niter,nsph*3)
      real*8 :: E0, E_plus, err, eeps, h
      integer :: iter, icomp, ksph, nsph_save, j, ncav_save, nbasis_save, iprint_save
!
      character(len=10) :: x1,x2
      character(len=30) :: fname
      integer :: fp, np
      logical :: iquiet_save
!
!---------------------------------------------------------------------------------------
!
!     save flags
      iquiet_save = iquiet ; iprint_save = iprint
!      
!     activate quiet mode
      iquiet = .true. ; iprint = 0
!
!     store esolv
      E0 = esolv
!
!     compute number of points in switch region
      call npoints_switch_region( np )
!
!     save initial DS
      nsph_save = nsph
      ncav_save = ncav
      nbasis_save = nbasis
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
!     set initial increment
      eeps=0.001d0
!      
!     initialize
      rwork = zero ; rrate = zero ; ework = zero ; fwork = zero ; eework = zero
!
!     loop over increments
      do iter = 1,niter
!        
!       loop over d / dr_k
        do ksph = 1,nsph_save
!
!         loop over components of d / dr_k
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph)*(1.d0+eeps) ; h = eeps*x_save(ksph)
            case(2) ; y(ksph) = y_save(ksph)*(1.d0+eeps) ; h = eeps*y_save(ksph)
            case(3) ; z(ksph) = z_save(ksph)*(1.d0+eeps) ; h = eeps*z_save(ksph)
            endselect
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!            
!           potential Phi and Psi vector
            call mkrhs( nsph_save, charge, x, y, z, ncav_save, ccav, phi, nbasis_save, psi )
!
!           solve PCM equations       
!!!            g(:,:)=zero ; sigma(:,:)=zero
!!!            call wghpot( phi, g )
!!!            call iefpcm( phi, g, psi, sigma, phi_eps, E_plus )
            call pcm(   .false.,  .true., .true., phi,      xx, phi_eps )
            call cosmo( .false., .false.,          xx, phi_eps, psi, sigma, E_plus )


!
!           compute energy
!!!            E_plus = 0.5d0 * (eps-1.d0)/eps * sprod( nbasis*nsph, sigma, psi )
!
!           account for a potentially null shift
            if ( abs(h) .gt. 1.E-12 ) then
!
!             account for a potentially null component of the force
              if ( abs( f(icomp,ksph) ).gt.1.E-12 ) then
!                        
!               compute relative error
!!!                err = abs( (E_plus - E0) / eeps + f(ksph,icomp) ) / abs( f(ksph,icomp) )
                err = abs( (E_plus - E0) / h + f(icomp,ksph) ) / abs( f(icomp,ksph) )
!
!               store
                rwork(iter,(ksph-1)*3+icomp) = err
                hwork(iter,(ksph-1)*3+icomp) = h
                ework(iter,(ksph-1)*3+icomp) = ( E_plus - E0 ) / h
                eework(iter,(ksph-1)*3+icomp) = E_plus
                fwork(iter,(ksph-1)*3+icomp) = f(icomp,ksph)

                
!!!                if ( abs( ( E_plus - E0 ) / h ) .gt. 1.E+03) then
!!!                  write(*,*) 'E_plus = ',E_plus
!!!                  write(*,*) 'E0     = ',E0
!!!                  stop
!!!                endif

!               compute rate
                if ( iter.gt.1 ) then 
                  rrate(iter,(ksph-1)*3+icomp) =  log( rwork(iter-1,(ksph-1)*3+icomp) / &
                                                       rwork(iter  ,(ksph-1)*3+icomp)   ) / &
                                                  log( hwork(iter  ,(ksph-1)*3+icomp) / &
                                                       hwork(iter-1,(ksph-1)*3+icomp)   )
                endif
              endif
            endif
!
          enddo
        enddo
!
        eeps = eeps / 2.d0
!
      enddo
!
!     print header
      write(*,*)'-----------------------------------------------------------------------------------------------'
      if ( np.gt.0 ) then 
        write(*,7324) lmax,ngrid
      else
        write(*,7325) lmax,ngrid
      endif
 7324 format(' lmax = ',i2,' , ngrid = ',i4,' - Integration points PRESENT in the switch region')     
 7325 format(' lmax = ',i2,' , ngrid = ',i4,' - Integration points NOT PRESENT in the switch region')     
      write(*,*)''
!      
!     print numerical derivative of energy
      write(*,8234) E0
 8234 format(' Energy (E_ref = ',e12.5,' ) : ')
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' E    r_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( eework(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''

      write(*,*)'Numerical derivative of energy : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( ework(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!      
!     print analytical force
      write(*,*)'Analytical force : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' force  _'i2','i1' : ',300(e12.5,2x))") j,icomp, ( fwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!     print relative error
      write(*,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!      
!     print rate of convergence
      write(*,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dE / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!
!     print to file
!     -------------
!
      write(x1,'(I2.2)') lmax
      write(x2,'(I4.4)') ngrid
      fname = 'lmax' // trim(x1) // '_ngrid' // trim(x2)
      fp    = 17
!
      open( unit=fp, file=fname, form='formatted', access='sequential', status='unknown')
!      
!     print header
      if ( np.gt.0 ) then 
        write(fp,7324) lmax,ngrid
      else
        write(fp,7325) lmax,ngrid
      endif
      write(fp,*)''
!
!
!     print numerical derivative of energy
      write(fp,*)'Numerical derivative of energy : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( ework(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!      
!     print analytical force
      write(fp,*)'Analytical force : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' force  _'i2','i1' : ',300(e12.5,2x))") j,icomp, ( fwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!
!     print relative error
      write(fp,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(fp,*) ''
!      
!     print rate of convergence
      write(fp,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(fp,"(' dE / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
!
      close(fp)
!      
!     restore DS
      call memfree
      call ddinit( nsph_save, x_save, y_save, z_save, r_save )
!
!     restore flags
      iquiet = iquiet_save ; iprint = iprint_save
!
!
endsubroutine check_forcesPCM
!-----------------------------------------------------------------------------------------
!
!
!
!
!-----------------------------------------------------------------------------------------
! Purpose : check derivatives of PCM matrix.
!-----------------------------------------------------------------------------------------
!
subroutine check_der_A()
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, zero, one, two, csph, rsph, memfree, &
                           ddinit, do_diag, eps
!
      implicit none
      real*8 :: f(nbasis,nsph)
      real*8 :: e(nbasis,nsph)
      real*8 :: A(nbasis*nsph,nbasis*nsph)
      real*8 :: dA(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: A_plus(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: xlm(nbasis),xx(ngrid),vplm(nbasis),vcos(lmax+1),vsin(lmax+1),basloc(nbasis)
      integer :: isph,jsph,i,j,ibeg,iend,nsph_save,icomp,ksph,iter,n
      real*8 :: eeps,err,rnorm
      integer, parameter :: niter = 6
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), s3(3)
      real*8 :: x(nsph), y(nsph), z(nsph), iwork(nsph*3,2), rwork(niter,nsph*3)
      real*8 :: rrate(niter,nsph*3)
      real*8 :: h,hwork(niter,nsph*3),rnorm0(nsph,nsph,nsph*3,niter), &
                                      err0(  nsph,nsph,nsph*3,niter)
      logical :: do_diag_save
!
!-----------------------------------------------------------------------------------------
!
!     save
      do_diag_save = do_diag
      do_diag = .true. 
!
!     build analytical derivative
!     ---------------------------
!
      dA(:,:,:,:) = zero
!
      do ksph = 1,nsph 
        do isph = 1,nsph
          do i = 1,nbasis
            do jsph = 1,nsph
              do j = 1,nbasis
!
!               standard basis vectors e_i, f_j
                e(:,:   )=zero
                f(:,:   )=zero
                e(i,isph)=one
                f(j,jsph)=one
!
!               compute < e, dA/dr_k f > 
                s3=zero
                call service_routine1_new( e(:,:), f(:,:), ksph, s3 )
!
!               store
                do icomp = 1,3
!
                  dA( (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp ) = s3(icomp)
!            
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
!
!     construct A
!     -----------
!
      A(:,:) = zero
!      
      do isph = 1,nsph
!
        do jsph = 1,nsph
          do j = 1,nbasis
!
!           standard basis vector e_j
            e(:,   :) = zero
            e(j,jsph) = one
!
!           compute A  _i e_j 
            ibeg = (isph-1)*nbasis+1
            iend = (isph-1)*nbasis+nbasis
!
!           recall to account for do_diag !!!
            call mkrvec( isph, eps, e, A( ibeg:iend,(jsph-1)*nbasis+j), xlm, xx, basloc, vplm, vcos, vsin )
!        
          enddo
        enddo
      enddo
!
!!!      do isph = 1,nsph
!!!        do i = 1,nbasis
!!!          write(*,"(4x,300(e12.5,2x))") ( A( (isph-1)*nbasis+i , j ) , j=1,nbasis*nsph )
!!!        enddo
!!!      enddo
!!!      write(*,*)''
!!!      read(*,*)
!
!
!     build FD approximation
!     ----------------------
!
!     save initial DS
      nsph_save = nsph
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
      rwork = zero ; hwork = zero ; rrate = zero ; err0 = zero ; rnorm0 = zero
!     set initial increment
      eeps=0.001d0
!
!     loop over increments
      do iter = 1,niter
!        
!       compute A^+
        do ksph = 1,nsph_save
!
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph)*(1.d0+eeps) ; h = eeps*x_save(ksph)
            case(2) ; y(ksph) = y_save(ksph)*(1.d0+eeps) ; h = eeps*y_save(ksph)
            case(3) ; z(ksph) = z_save(ksph)*(1.d0+eeps) ; h = eeps*z_save(ksph)
            endselect
!
!           if null increment, skip case
            if ( abs(h).lt.1.E-12 )  cycle
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!
!           initialize
            A_plus(:,:,ksph,icomp)=zero
!
!           build A^+
            do isph = 1,nsph
              do jsph = 1,nsph
                do j = 1,nbasis
!
                  e(:,:   )=zero
                  e(j,jsph)=one
!
!                 compute A^+_i e_j
                  ibeg = (isph-1)*nbasis+1
                  iend = (isph-1)*nbasis+nbasis
!
!                 recall to account for do_diag !!!
                  call mkrvec( isph, eps, e, A_plus( ibeg:iend,(jsph-1)*nbasis+j,ksph,icomp ), xlm, xx, basloc, vplm, vcos, vsin )
!
!!!!                 accumulate error
!!!                  do i = 1,nbasis
!!!                    err = err + ( ( A_plus((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) -        &
!!!                                    A(     (isph-1)*nbasis+i,(jsph-1)*nbasis+j           ) ) / h -  &
!!!                                    dA(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)           )**2
!!!!
!!!                  enddo
!
                enddo
              enddo
            enddo
!
!           accumulate error
            err=zero ; rnorm=zero 
!            
            do isph =1,nsph
              do jsph =1,nsph
!                
                do i = 1,nbasis
                  do j = 1,nbasis
!                  
!                  
!                   sphere-to-sphere relative error 
                    err0(isph,jsph,(ksph-1)*3+icomp,iter) = &
                    err0(isph,jsph,(ksph-1)*3+icomp,iter) + &
!                    
                    ( ( A_plus((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) -        &
                        A(     (isph-1)*nbasis+i,(jsph-1)*nbasis+j           ) ) / h -  &
                        dA(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)           )**2
!                                    
!                                    
!                   sphere-to-sphere norm
                    rnorm0(isph,jsph,(ksph-1)*3+icomp,iter) = &
                    rnorm0(isph,jsph,(ksph-1)*3+icomp,iter) + &
!                    
                    dA((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)**2
!
                  enddo
                enddo
!                
                err   = err   + err0(  isph,jsph,(ksph-1)*3+icomp,iter)
                rnorm = rnorm + rnorm0(isph,jsph,(ksph-1)*3+icomp,iter)
!                
              enddo
            enddo
!
!           take square root
            err = sqrt(err)
            rnorm = sqrt(rnorm)

            do isph = 1,nsph
              do jsph = 1,nsph

                err0(isph,jsph,(ksph-1)*3+icomp,iter) = &
                sqrt( err0(isph,jsph,(ksph-1)*3+icomp,iter)/rnorm )
!!!                sqrt( err0(isph,jsph,(ksph-1)*3+icomp,iter)/rnorm / &
!!!                      rnorm0(isph,jsph,(ksph-1)*3+icomp,iter)/rnorm )

              enddo
            enddo
!
!           store
            rwork(iter,(ksph-1)*3+icomp) = err/rnorm
            hwork(iter,(ksph-1)*3+icomp) = h
!
!!!            do isph = 1,nsph
!!!              write(*,"(' err('i2',:) = ', 300(e12.5,2x))") isph, ( sqrt(err0(isph,jsph))/rnorm , jsph=1,nsph )
!!!            enddo
!!!            read(*,*)
!
!           compute rate
            if ( iter.gt.1 ) then 
              rrate(iter,(ksph-1)*3+icomp) =  log( rwork(iter-1,(ksph-1)*3+icomp) / &
                                                   rwork(iter  ,(ksph-1)*3+icomp)   ) / &
                                              log( hwork(iter  ,(ksph-1)*3+icomp) / &
                                                   hwork(iter-1,(ksph-1)*3+icomp)   )
            endif
! 
          enddo
        enddo
!
!       update
        eeps = eeps / 2.d0
!        
      enddo
!
!     extensive printing
      do ksph=1,nsph
        do icomp=1,3
!
          write(*,*)'---------------------------------------------------------------'
          write(*,"(' dA / dr_'i2','i1' rel. error : ')") ksph,icomp
          write(*,*)''
!
          do iter=1,niter
            do isph = 1,nsph

              write(*,"(' [',300(e12.5,2x),' ]')") ( err0(isph,jsph,(ksph-1)*3+icomp,iter) , jsph=1,nsph )

            enddo
            write(*,*)''
          enddo
        enddo
      enddo
      read(*,*)
!
!     print relative error
      write(*,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dA / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!      
!     print rate of convergence
      write(*,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dA / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
      read(*,*)
!
!     restore
      do_diag = do_diag_save
!
!
endsubroutine check_der_A
!----------------------------------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------------------------------
! Purpose : check derivatives of indicator function U_j^n.
!----------------------------------------------------------------------------------------------
!
subroutine check_der_U()
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, zero, one, two, csph, rsph, memfree, &
                           ddinit, ui, du
!
      implicit none
      integer :: isph,jsph,i,j,ibeg,iend,nsph_save,icomp,ksph,iter,n
      real*8 :: eeps,rnorm
      integer, parameter :: niter = 6
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph)
      real*8 :: x(nsph), y(nsph), z(nsph), rwork(niter,nsph*3,nsph)
      real*8 :: rrate(niter,nsph*3,nsph)
      real*8 :: h,hwork(niter,nsph*3,nsph)
      real*8 :: u0(ngrid,nsph),du0(3,nsph,ngrid,nsph),error
!
!-----------------------------------------------------------------------------------------
!
!     store
      u0(:,:) = ui(:,:)
      du0(:,:,:,:) = du(:,:,:,:)
!
!     build FD approximation
!     ----------------------
!
!     save initial DS
      nsph_save = nsph
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
      rwork = zero ; hwork = zero ; rrate = zero 
!     set initial increment
      eeps=0.001d0
!
!     loop over increments
      do iter = 1,niter
!        
!       compute U+
        do ksph = 1,nsph_save
!
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph)*(1.d0+eeps) ; h = eeps*x_save(ksph)
            case(2) ; y(ksph) = y_save(ksph)*(1.d0+eeps) ; h = eeps*y_save(ksph)
            case(3) ; z(ksph) = z_save(ksph)*(1.d0+eeps) ; h = eeps*z_save(ksph)
            endselect
!
!           if null increment, skip case
            if ( abs(h).lt.1.E-12 )  cycle
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!
!           loop over spheres
            do isph = 1,nsph
!
!             initialize
              error = zero ; rnorm = zero
!
!             loop over integration points
              do i = 1,ngrid
!
!               accumulate
                error = error + ( ( ui(i,isph) - u0(i,isph) )/h - du0(icomp,ksph,i,isph) )**2
                rnorm = rnorm + (                                 du0(icomp,ksph,i,isph) )**2
!
              enddo
!
!             take square root
              error = sqrt(error)
              rnorm = sqrt(rnorm)

!             store
              if ( rnorm.ge.1.E-10 ) then
              rwork(iter,(ksph-1)*3+icomp,isph) = error/rnorm
              endif
              hwork(iter,(ksph-1)*3+icomp,isph) = h
!              
!             compute rate
              if ( iter.gt.1 ) then 
                if ( rwork(iter,(ksph-1)*3+icomp,isph) .ge. 1.E-10 ) then
              rrate(iter,(ksph-1)*3+icomp,isph) =  log( rwork(iter-1,(ksph-1)*3+icomp,isph) / &
                                                        rwork(iter  ,(ksph-1)*3+icomp,isph)   ) / &
                                                   log( hwork(iter  ,(ksph-1)*3+icomp,isph) / &
                                                        hwork(iter-1,(ksph-1)*3+icomp,isph)   )
              endif
              endif
!
            enddo
          enddo
        enddo
!
!       update
        eeps = eeps / 2.d0
!        
      enddo
!
!     extensive printing
      write(*,*) 'Relative error :'
      do isph=1,nsph
        do ksph=1,nsph
          do icomp=1,3
!
            write(*,"(' dU_'i2' / dr_'i2','i1'   ',300(e12.5,2x))") isph,ksph,icomp, &
            ( rwork(iter,(ksph-1)*3+icomp,isph) , iter=1,niter )
!
          enddo
        enddo
      enddo
      write(*,*) ''
!      
!     print rate of convergence
      write(*,*) 'Rate of convergence :'
      do isph=1,nsph
        do ksph=1,nsph
          do icomp=1,3
!
            write(*,"(' dU_'i2' / dr_'i2','i1'   ',300(e12.5,2x))") isph,ksph,icomp, &
            ( rrate(iter,(ksph-1)*3+icomp,isph) , iter=1,niter )
!
          enddo
        enddo
      enddo
      write(*,*) ''
!
!
endsubroutine check_der_U
!-----------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------
! Purpose : check derivatives of COSMO matrix.
!-----------------------------------------------------------------------------------------
!
subroutine check_der_L()
!
      use ddcosmo , only : nbasis, nsph, ngrid, lmax, zero, one, two, csph, rsph, memfree, &
                           ddinit, fdoka, fdokb, basis, do_diag
!
      implicit none
      real*8 :: f(nbasis,nsph)
      real*8 :: e(nbasis,nsph),Ae(nbasis,nsph)
      real*8 :: A(nbasis*nsph,nbasis*nsph)
      real*8 :: dA(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: A_plus(nbasis*nsph,nbasis*nsph,nsph,3)
      real*8 :: xlm(nbasis),xx(ngrid),vplm(nbasis),vcos(lmax+1),vsin(lmax+1), &
                basloc(nbasis),dbsloc(3,nbasis),xi(ngrid,nsph)
      integer :: isph,jsph,i,j,ibeg,iend,nsph_save,icomp,ksph,iter,n,k
      real*8 :: eeps,err,rnorm
      integer, parameter :: niter = 6
      real*8 :: x_save(nsph), y_save(nsph), z_save(nsph), r_save(nsph), s3(3)
      real*8 :: x(nsph), y(nsph), z(nsph), rwork(niter,nsph*3)
      real*8 :: rrate(niter,nsph*3)
      real*8 :: h,hwork(niter,nsph*3),rnorm0(nsph,nsph,nsph*3,niter), &
                                      err0(  nsph,nsph,nsph*3,niter)
      logical :: do_diag_save
!
!-----------------------------------------------------------------------------------------
!
      do_diag_save = do_diag
      do_diag = .true.
!      
!     build analytical derivative dL
!     ------------------------------
!
      dA(:,:,:,:) = zero
!
      do ksph = 1,nsph 
        do isph = 1,nsph
          do i = 1,nbasis
!          
!           standard basis vector e_i
            e(:,:   )=zero
            e(i,isph)=one
!
!           expand e_i in spherical harmonics
            xi(:,:) = zero
            do k = 1,nsph
              do n = 1,ngrid
!              
                xi(n,k) = dot_product( e(:,k), basis(:,n) )
!                
              enddo
            enddo
!                
            do jsph = 1,nsph
              do j = 1,nbasis
!
!               standard basis vector f_j
                f(:,:   )=zero
                f(j,jsph)=one
!                
!               compute < e, dL/dr_k f > 
                s3=zero
!                
                call fdoka( ksph, f, xi(:,ksph), basloc, dbsloc, vplm, vcos, vsin, s3 ) 
                call fdokb( ksph, f, xi,         basloc, dbsloc, vplm, vcos, vsin, s3 ) 
!
!               store
                do icomp = 1,3
!
                  dA( (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp ) = s3(icomp)
!            
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
!
!     construct L
!     -----------
!
      A(:,:) = zero
!      
      do isph = 1,nsph
!
        do jsph = 1,nsph
          do j = 1,nbasis
!
!           standard basis vector e_j
            e(:,:   ) = zero
            e(j,jsph) = one
!
!           compute L e_j
!           [ recall to account for do_diag !!! "j" is placeholder !!! ] 
            call Lx( j, e, Ae )
!
!           store Le
            do k = 1,nsph
!            
              ibeg = (k-1)*nbasis+1
              iend = (k-1)*nbasis+nbasis
!              
              A( ibeg:iend,(jsph-1)*nbasis+j) = Ae(:,k)
!              
            enddo
!        
          enddo
        enddo
      enddo
!
!
!     build FD approximation
!     ----------------------
!
!     save initial DS
      nsph_save = nsph
      x_save = csph(1,:)
      y_save = csph(2,:)
      z_save = csph(3,:)
      r_save = rsph(  :)
!
      rwork = zero ; hwork = zero ; rrate = zero ; err0 = zero ; rnorm0 = zero
!     set initial increment
      eeps=0.001d0
!
!     loop over increments
      do iter = 1,niter
!        
!       compute L^+
        do ksph = 1,nsph_save
!
          do icomp = 1,3
!
!           deallocate DS      
            call memfree
!
!           perturb     
            x = x_save
            y = y_save
            z = z_save
            select case(icomp)
            case(1) ; x(ksph) = x_save(ksph)*(1.d0+eeps) ; h = eeps*x_save(ksph)
            case(2) ; y(ksph) = y_save(ksph)*(1.d0+eeps) ; h = eeps*y_save(ksph)
            case(3) ; z(ksph) = z_save(ksph)*(1.d0+eeps) ; h = eeps*z_save(ksph)
            endselect
!
!           if null increment, skip case
            if ( abs(h).lt.1.E-12 )  cycle
!
!           allocate new DS      
            call ddinit( nsph_save, x, y, z, r_save )
!
!           initialize
            A_plus(:,:,ksph,icomp)=zero
!
!           build L^+
            do isph = 1,nsph
              do jsph = 1,nsph
                do j = 1,nbasis
!
                  e(:,:   )=zero
                  e(j,jsph)=one
!
!                 compute L e_j
!                 [ recall to account for do_diag !!! "j" is placeholder !!! ] 
                  call Lx( j, e, Ae )
!
!                 store Le
                  do k = 1,nsph
!
                    ibeg = (k-1)*nbasis+1
                    iend = (k-1)*nbasis+nbasis
!                    
                    A_plus( ibeg:iend,(jsph-1)*nbasis+j,ksph,icomp ) = Ae(:,k) 
!                    
                  enddo
!
!!!!                 accumulate error
!!!                  do i = 1,nbasis
!!!                    err = err + ( ( A_plus((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) -        &
!!!                                    A(     (isph-1)*nbasis+i,(jsph-1)*nbasis+j           ) ) / h -  &
!!!                                    dA(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)           )**2
!!!!
!!!                  enddo
!
                enddo
              enddo
            enddo
!
!           accumulate error
            err=zero ; rnorm=zero 
!            
            do isph =1,nsph
              do jsph =1,nsph
!                
                do i = 1,nbasis
                  do j = 1,nbasis
!                  
!                  
!                   sphere-to-sphere relative error 
                    err0(isph,jsph,(ksph-1)*3+icomp,iter) = &
                    err0(isph,jsph,(ksph-1)*3+icomp,iter) + &
!                    
                    ( ( A_plus((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp) -        &
                        A(     (isph-1)*nbasis+i,(jsph-1)*nbasis+j           ) ) / h -  &
                        dA(    (isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)           )**2
!                                    
!                                    
!                   sphere-to-sphere norm
                    rnorm0(isph,jsph,(ksph-1)*3+icomp,iter) = &
                    rnorm0(isph,jsph,(ksph-1)*3+icomp,iter) + &
!                    
                    dA((isph-1)*nbasis+i,(jsph-1)*nbasis+j,ksph,icomp)**2
!
                  enddo
                enddo
!                
                err   = err   + err0(  isph,jsph,(ksph-1)*3+icomp,iter)
                rnorm = rnorm + rnorm0(isph,jsph,(ksph-1)*3+icomp,iter)
!                
              enddo
            enddo
!
!           take square root
            err = sqrt(err)
            rnorm = sqrt(rnorm)

            do isph = 1,nsph
              do jsph = 1,nsph

                err0(isph,jsph,(ksph-1)*3+icomp,iter) = &
                sqrt( err0(isph,jsph,(ksph-1)*3+icomp,iter)/rnorm )
!!!                sqrt( err0(isph,jsph,(ksph-1)*3+icomp,iter)/rnorm / &
!!!                      rnorm0(isph,jsph,(ksph-1)*3+icomp,iter)/rnorm )

              enddo
            enddo
!
!           store
            rwork(iter,(ksph-1)*3+icomp) = err/rnorm
            hwork(iter,(ksph-1)*3+icomp) = h
!
!!!            do isph = 1,nsph
!!!              write(*,"(' err('i2',:) = ', 300(e12.5,2x))") isph, ( sqrt(err0(isph,jsph))/rnorm , jsph=1,nsph )
!!!            enddo
!!!            read(*,*)
!
!           compute rate
            if ( iter.gt.1 ) then 
              rrate(iter,(ksph-1)*3+icomp) =  log( rwork(iter-1,(ksph-1)*3+icomp) / &
                                                   rwork(iter  ,(ksph-1)*3+icomp)   ) / &
                                              log( hwork(iter  ,(ksph-1)*3+icomp) / &
                                                   hwork(iter-1,(ksph-1)*3+icomp)   )
            endif
! 
          enddo
        enddo
!
!       update
        eeps = eeps / 2.d0
!        
      enddo
!
!     extensive printing
      do ksph=1,nsph
        do icomp=1,3
!
          write(*,*)'---------------------------------------------------------------'
          write(*,"(' dL / dr_'i2','i1' rel. error : ')") ksph,icomp
          write(*,*)''
!
          do iter=1,niter
            do isph = 1,nsph

              write(*,"(' [',300(e12.5,2x),' ]')") ( err0(isph,jsph,(ksph-1)*3+icomp,iter) , jsph=1,nsph )

            enddo
            write(*,*)''
          enddo
        enddo
      enddo
      read(*,*)
!
!     print relative error
      write(*,*)'Relative error : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dL / dr_'i2','i1' : ',300(e12.5,2x))") j,icomp, ( rwork(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
!      
!     print rate of convergence
      write(*,*)'Rate of convergence : '
      do j = 1,nsph
        do icomp = 1,3
!
          write(*,"(' dL / dr_'i2','i1' : ',300(f12.3,2x))") j,icomp, ( rrate(iter,(j-1)*3+icomp) , iter=1,niter )
!        
        enddo
      enddo
      write(*,*) ''
      read(*,*)
!
!
      do_diag = do_diag_save
!
!
endsubroutine check_der_L
!
!
!-----------------------------------------------------------------------------------------
subroutine npoints_switch_region( np )
!
      use ddcosmo , only : ngrid, nsph, csph, rsph, grid, inl, nl, eta, se
!
      implicit none
!
      integer, intent(out) :: np
!      
      real*8  :: v_nij(3), t_nij, upper, lower
      integer :: n, isph, jsph, j
!
!-----------------------------------------------------------------------------------------
!
!     intialize point counter
      np = 0

!     loop over integration points
      do n = 1,ngrid
!
!       loop over spheres
        do isph = 1,nsph
!
!         loop over neighbors of i-sphere
          do j = inl(isph),inl(isph+1)-1
!
!           neighbor is j-sphere            
            jsph = nl(j)
!
!           t
            v_nij = csph(:,isph) + rsph(isph)*grid(:,n) - csph(:,jsph)
            t_nij = sqrt( dot_product( v_nij,v_nij ) ) / rsph(jsph)
!
!           upper and lower bounds of switch region
            upper = 1.d0 + (se + 1.d0)*eta / 2.d0
            lower = upper - eta
!
!           increment counter
            if ( ( t_nij.ge.lower ) .and. ( t_nij.le.upper ) )  np = np + 1
!
          enddo
        enddo
      enddo
!
!
endsubroutine npoints_switch_region
!-----------------------------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------------------------
!
subroutine test_fmm()
!
      use ddcosmo , only : nsph, eps, nbasis, ngrid, lmax, do_diag
!
      implicit none
      real*8  :: vlm(nbasis,nsph), dvlm(nbasis), x(ngrid), xlm(nbasis), basloc(nbasis), & 
                 vplm(nbasis), vcos(lmax+1), vsin(lmax+1), dvlm_fmm(nbasis), &
                 error, rnorm
      integer :: isph,i
!      
!-----------------------------------------------------------------------------------------

      write(*,*)'FMM is not working yet !'
      return
!
!     initialize random number generator
      call random_seed
!      
!     build a random vector vlm
      do isph = 1,nsph
        do i = 1,nbasis
!
          call random_number( vlm(i,isph) )
!
        enddo
      enddo
!      
!     loop over spheres
      do isph = 1,nsph
!      
        do_diag = .true.
        call mkrvec(     isph, eps, vlm, dvlm,     xlm, x, basloc, vplm, vcos, vsin )
        call mkrvec_fmm( isph, eps, vlm, dvlm_fmm, xlm, x, basloc, vplm, vcos, vsin )
!      
!       initialize
        error = 0.d0 ; rnorm = 0.d0
!
        do i = 1,nbasis
!
!         accumulate
          error = error + ( dvlm(i) - dvlm_fmm(i) )**2
          rnorm = rnorm + ( dvlm(i)               )**2
!          
        enddo
!        
      enddo
!        
      write(*,1000) sqrt(error/rnorm)
 1000 format(' relative error = ',e12.5)
!
!
endsubroutine test_fmm
!-----------------------------------------------------------------------------------------
!
!
!
!
!
!-----------------------------------------------------------------------------------------
subroutine check_old_new_pcm()
!
      use ddcosmo , only : ngrid, nsph, wghpot, ncav, nbasis, memfree, &
                           read_x, read_y, read_z, read_r, read_q, ddinit, &
                           read_molecule_file, ccav, tokcal, iquiet, iprint, &
                           use_fmm
!
      implicit none
      real*8 :: g(ngrid,nsph), phi(ncav), psi(nbasis,nsph), sigma(nbasis,nsph), &
                phi_eps(nbasis,nsph), esolv, xx(1), charge(nsph)
      integer :: iprint_save
      logical :: iquiet_save, use_fmm_save
!
!-----------------------------------------------------------------------------------------
!
      use_fmm_save = use_fmm
!
      iquiet_save = iquiet ; iprint_save = iprint
      iquiet = .true. ; iprint = 0
!      
!      
!     LEGACY, FMM
!     ===========

!     free the memory
      call memfree()
!
      use_fmm = .true.
!      
!     recreate nsph, read_x, read_y, read_z, read_r, read_q
      call read_molecule_file()
!
!     redirect charge [ THIS WILL NEED TO BE CLEANED UP ]
      charge = read_q
!
!     initialize datastructure
      call ddinit( nsph, read_x, read_y, read_z, read_r )
!      
!     solute electrostatic potential phi, and psi vector
      call mkrhs( nsph, charge, read_x, read_y, read_z, ncav, ccav, phi, nbasis, psi )
!
!     pcm
      call wghpot( phi, g )
      call iefpcm( phi, g, psi, sigma, phi_eps, esolv )
      write (6,'(1x,a,f14.6)') '[OLD,FMM] electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
!      
!     LEGACY, NO FMM
!     ==============

!     free the memory
      call memfree()
!
      use_fmm = .false.
!      
!     recreate nsph, read_x, read_y, read_z, read_r, read_q
      call read_molecule_file()
!
!     redirect charge [ THIS WILL NEED TO BE CLEANED UP ]
      charge = read_q
!
!     initialize datastructure
      call ddinit( nsph, read_x, read_y, read_z, read_r )
!      
!     solute electrostatic potential phi, and psi vector
      call mkrhs( nsph, charge, read_x, read_y, read_z, ncav, ccav, phi, nbasis, psi )
!
!     pcm
      call wghpot( phi, g )
      call iefpcm( phi, g, psi, sigma, phi_eps, esolv )
      write (6,'(1x,a,f14.6)') '[OLD    ] electrostatic solvation energy (kcal/mol):', esolv*tokcal
!
!
!     NEW
!     ===

!     free the memory
      call memfree()
!      
!     recreate nsph, read_x, read_y, read_z, read_r, read_q
      call read_molecule_file()
!
!     redirect charge [ THIS WILL NEED TO BE CLEANED UP ]
      charge = read_q
!
!     initialize datastructure
      call ddinit( nsph, read_x, read_y, read_z, read_r )
!      
!     solute electrostatic potential phi, and psi vector
      call mkrhs( nsph, charge, read_x, read_y, read_z, ncav, ccav, phi, nbasis, psi )
!
      call pcm(   .false.,  .true., .true., phi,      xx, phi_eps )
      call cosmo( .false., .false.,          xx, phi_eps, psi, sigma, esolv )
      write (6,'(1x,a,f14.6)') '[NEW    ] electrostatic solvation energy (kcal/mol):', esolv*tokcal
!      
!
!     clean up
!     ========
!
!     free the memory
      call memfree()
!      
!     recreate nsph, read_x, read_y, read_z, read_r, read_q
      call read_molecule_file()
!      
!     initialize datastructure
      call ddinit( nsph, read_x, read_y, read_z, read_r )
!
      iquiet = iquiet_save ; iprint = iprint_save ; use_fmm = use_fmm_save
!
!
endsubroutine
