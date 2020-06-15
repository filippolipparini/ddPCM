subroutine incore
  use ddcosmo
  implicit none
!
  integer :: isph, memreq, blkmem, istat, lnl, ij, jsph, its, l, ind, mm, m
  real*8  :: vij(3), sij(3), vvij, tij, xij, oij, tt, ft
  real*8,  allocatable :: basloc(:), vplm(:), vcos(:), vsin(:), scr(:,:), blk(:,:)
  integer, allocatable :: iscr(:)
!
!
! compute the required amount of memory to store the matrix incore:
!
  blkmem = nylm*nylm
  lnl    = inl(nsph+1) - 1
  memreq = lnl*blkmem
!
  allocate (lmat(nylm,nylm,lnl),stat=istat)
  allocate (sig_gather(nylm,maxnl))
  havemem = istat.eq.0
!
  if (havemem) then
    write(6,*) ' memory used for the matrix:', float(8*memreq)/1024.d0/1024.d0
!
! there is enough memory to build the l matrix incore. 
! do it!
!
!   lmat = zero
    allocate (basloc(nylm),vplm(nylm),vcos(lmax+1),vsin(lmax+1))                                                                      
    allocate (scr(nylm,ngrid),blk(nylm,nylm))
!
!   create an index of blocks for the transpose matrix:
!
    allocate (iltrn(lnl),iscr(nsph))
!
    iscr = 0
    do isph = 1, nsph
      do ij = inl(isph), inl(isph+1) - 1
        jsph = nl(ij)
        iltrn(inl(jsph)+iscr(jsph)) = ij
        iscr(jsph) = iscr(jsph) + 1
      end do
    end do
!
    deallocate (iscr)
!
!$omp parallel do default(shared) private(isph,ij,jsph,its,vij,vvij,tij,sij,xij,oij,tt,l,ind,ft,mm,m) &
!$omp private(scr,basloc,vplm,vcos,vsin,blk) schedule(guided)
    do isph = 1, nsph
      do ij = inl(isph), inl(isph+1) - 1
        jsph = nl(ij) 
        scr  = zero
        do its = 1, ngrid
          if (ui(its,isph).lt.one) then
            vij  = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
            vvij = sqrt(vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3))                                                                    
            tij  = vvij/rsph(jsph)                                                                                                        
!
            if (tij.lt.one) then
!
              sij = vij/vvij
              xij = fsw(tij, se, eta)
              if(fi(its,isph).gt.one) then
                oij = xij/fi(its,isph)
              else                                                                                                                        
                oij = xij                                                                                                                 
              end if                                                                                                                      
              call ylmbas(sij,basloc,vplm,vcos,vsin)
!
!             integrate:
!
              tt = oij*w(its)
              do l = 0, lmax
                ind = l*l + l + 1
                ft = tt/facl(ind)
                do mm = 0, 2*l
                  m = mm - l
                  scr(ind+m,its) = scr(ind+m,its) + ft*basloc(ind+m)
                end do
                tt = tt*tij
              end do
            end if
          end if
        end do
!       call dgemm('N','T',nylm,nylm,ngrid,one,basis,nylm,scr,nylm,zero,lmat(1,1,ij),nylm)
        call dgemm('N','T',nylm,nylm,ngrid,one,basis,nylm,scr,nylm,zero,blk,nylm)
        lmat(:,:,ij) = blk
!       write(6,*) 'isph, ij, blk:', isph, ij
!       write(6,'(25f12.6)') lmat(:,:,ij)
      end do
    end do
!
    deallocate (basloc, vplm, vcos, vsin, scr)
  else
    return
  end if
!
  return
end subroutine incore
