!--------------------------------------------------------------------------------------
! Purpose : 
!
! 1) Collocate rhs of eq. (14) at LL integration points.
!
! Rhs :
!
!   - ( 1 - sum \omega^jk(s) ) Phi_j(s)   ,   s \in S^2
!            k                    
!
! For each domain \Omega_j, collocate at {s_n} : 
!
!   phi(n,j) = - ( 1 - sum \omega^jk(s_n) ) Phi_j(s_n) 
!                       k         
!
! Rename :
!
!   Phi_n^j = Phi_j(s_n)
!
!   W_n^jk = \omega^jk(s_n)
!
!   U_n^j = 1 - sum W_n^jk
!                k
! 
! Thus :
!
!   phi(n,j) = - U_n^j * Phi_n*j
!
! Remark : quantities phi and Phi should not be confused!
!
! 2) For each atom, compute
!
!   psi(i) = 2pi^1/2  q_i
!
!--------------------------------------------------------------------------------------
!
!
subroutine dorhs( doen, phi, psi, en )
!
      use  ddcosmo
!      
      implicit none
      logical,                        intent(in)    :: doen
      real*8, dimension(ngrid,nsph),  intent(inout) :: phi, en
      real*8, dimension(nbasis,nsph), intent(inout) :: psi
!
      integer :: isph, its, j
      real*8  :: v, ex, ey, ez
      real*8  :: dx, dy, dz, d2, d, d3, fe, fac
!
!--------------------------------------------------------------------------------------
!
!     compute scaling factor
      fac = sqrt(four*pi)
!
!     initialize
      phi = zero
      psi = zero
      if (doen) en = zero
!      
      !
      !$omp parallel do default(shared) private(isph,its,j,v,ex,ey,ez,dx,dy,dz,d2,d3,d,fe)
!
!     1st loop over atoms
      do isph = 1, nsph
!      
        psi(1,isph) = fac*q(isph)
!
!       loop over integration points
        do its = 1, ngrid
!
!         non-null contribution
          if (ui(its,isph).gt.zero) then
!                  
!           initialize electrostatic potential Phi
            v=zero
!            
            ex=zero ; ey=zero ; ez=zero
!
!           2nd loop over atoms
            do j = 1, nsph
!
!             d = | rr_j + r_j ss_n - rr_k |
              dx = csph(1,isph) + rsph(isph)*grid(1,its) - csph(1,j)
              dy = csph(2,isph) + rsph(isph)*grid(2,its) - csph(2,j)
              dz = csph(3,isph) + rsph(isph)*grid(3,its) - csph(3,j)
              d2 = dx*dx + dy*dy + dz*dz
              d  = sqrt(d2)
!              
!             d3 = d^3/2
              d3 = d*d2
              fe = q(j)/d3
!
!             accumulate for electrostatic potential Phi 
              v = v + q(j)/d
!              
              ex = ex + fe*dx
              ey = ey + fe*dy
              ez = ez + fe*dz
!
            end do
!            
!           assemble rhs entry
            phi(its,isph) = - ui(its,isph)*v
!            
      !     if (doen) en(its,isph) = - (ex*grid(1,its) + ey*grid(2,its) + ez*grid(3,its))
            if (doen) en(its,isph) = - ui(its,isph)*(ex*grid(1,its) + ey*grid(2,its) + ez*grid(3,its))
!            
          end if
        end do
      end do
!
!     printing
      if (iprint.ge.3) then
!
        call ptcart( 'potential', nsph, 0, phi )
        call ptcart( 'en'       , nsph, 0, en  )
!
      end if
      return
!
!
endsubroutine dorhs
