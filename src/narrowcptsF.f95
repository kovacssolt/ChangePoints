subroutine narrowcptsF(x,n,dec,dep,ilen,nint,nsum,par,thr,t,bou,bou2,cpts2) bind(c,name="narrowcptsF_")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      real(kind=c_double), parameter :: eps=1e-4
      integer(kind= c_int), intent(in)          :: n,dep,nsum,par
      integer(kind= c_int), intent(in)          :: nint(dep-1)
      real(kind=c_double),  intent(in)          :: x(n),dec,ilen(dep-1),thr
      integer(kind= c_int)                      :: i,j,k,l,m,ipm,prev,res1,low,up,j2,t,ord(nsum),t2(dep),t3,i2,i3
      real(kind=c_double)                       :: cu(n+1),ran(dep-1),temp(n-1),sl,sr,sf,prevf,res2,n2,help(nsum)
      integer(kind= c_int),  dimension (:), allocatable :: s(:,:)
      real(kind= c_double), intent(out)         :: bou(nsum)
      integer(kind= c_int), intent(out)         :: bou2(nsum,3),cpts2(n)
10000 j=2
      t2=0
      cu(1)=0
10001 do 10002 i=1,n
      cu(i+1)=cu(i)+x(i)
10002 continue
      bou2(1,1)=1
      bou2(1,2)=n
      sf=cu(n+1)
      if(par .lt. n) then
        prev=1+floor((n-1)/3.0_c_double)
        sl=cu(prev+1)
        sr=cu(n+1)-cu(prev+1)
        prevf=sqrt( sl**2 / (prev+1-bou2(1,1)) + sr**2 / (bou2(1,2)-prev) - sf**2 / (bou2(1,2)-bou2(1,1)+1) )
        call obsnarrowF(n, bou2(1,1), bou2(1,2), bou2(1,1), bou2(1,2), prev, cu, prevf,res1,res2)
        bou2(1,3)=res1
        bou(1)=res2
      else
10005 do 10006 l=1,n-1
        t=l+1
        sl=cu(t)
        sr=cu(n+1)-cu(t)
        temp(l)=sl**2/l + sr**2/(n-l)
10006 continue
        bou2(1,3)=maxloc(temp,1)
        bou(1)=sqrt(temp(bou2(1,3)) - sf**2/n)
      endif
      if(bou(1) .ge. thr) then
        t2(1)=1
      end if
10007 do 10014 i=1,dep-1
      ran(i) = (n-ilen(i))/(nint(i)-1)
      if(par .lt. ilen(i)+1) then
10008 do 10009 k=1,nint(i)
        bou2(j,1)=floor(1+(k-1)*ran(i))
        bou2(j,2)=ceiling(ilen(i)+(k-1)*ran(i)-1e-14)
        sf=cu(bou2(j,2)+1) - cu(bou2(j,1))
        prev=bou2(j,1)+floor((bou2(j,2)-bou2(j,1))/3.0_c_double)
        sl=cu(prev+1)-cu(bou2(j,1))
        sr=cu(bou2(j,2)+1)-cu(prev+1)
        prevf=sqrt( sl**2 / (prev+1-bou2(j,1)) + sr**2 / (bou2(j,2)-prev) - sf**2 / (bou2(j,2)-bou2(j,1)+1) )
        call obsnarrowF(n, bou2(j,1), bou2(j,2), bou2(j,1), bou2(j,2), prev, cu, prevf,res1,res2)
        bou2(j,3)=res1
        bou(j)=res2
        if(bou(j) .lt. thr) then
          !bou(j)=0
          help(j)=0
        else
          help(j)=bou(j)
          t2(i+1)=t2(i+1)+1
        end if
        j=j+1
10009 continue
      else
10010 do 10013 k=1,nint(i)
        bou2(j,1)=floor(1+(k-1)*ran(i))
        bou2(j,2)=ceiling(ilen(i)+(k-1)*ran(i)-1e-14)
        sf=cu(bou2(j,2)+1) - cu(bou2(j,1))
        m=bou2(j,2)-bou2(j,1)
10011 do 10012 l=1,m
        t = bou2(j,1)+l
        sl= cu(t)-cu(bou2(j,1))
        sr= cu(bou2(j,2)+1)-cu(t)
        temp(l)=sl**2/(t-bou2(j,1)) + sr**2/(bou2(j,2)-t+1)
10012 continue
        ipm=maxloc(temp(1:m),1)
        bou2(j,3)=ipm+bou2(j,1)-1
        bou(j)=sqrt( temp(ipm) - sf**2/(bou2(j,2)-bou2(j,1)+1) )
        if(bou(j) .lt. thr) then
          help(j)=0
        else
          help(j)=bou(j)
          t2(i+1)=t2(i+1)+1
        end if
        j=j+1
10013 continue
      end if
10014 continue
      cpts2=0
      t=0
      j2=nsum
10019 do 10029 i3=1,dep-1
      i2=dep-i3
      if(t2(i2+1) .ge. 1) then
10020 do 10028 k=1,nint(i2)
        if(help(j2-k) .ge. thr) then
          t=t+1
          cpts2(t)=bou2(j2-k,3)
          if(t .eq. 1) help(1)=0
          j=1
10021 do 10024 i=1,i2
          low=min(ceiling((cpts2(t)-1)/ran(i)+1),nint(i))
          if( cpts2(t) .lt. floor(1+(low-1)*ran(i))) then
            low=low-1
          else if( cpts2(t) .ge. floor(1+low*ran(i))) then
            low=min(low+1,nint(i))
          end if
          up=max(floor((cpts2(t)+1-ilen(i))/ran(i)+1),1)
          if( cpts2(t)+1 .gt. ceiling(ilen(i)+(up-1)*ran(i))) then
            up=up+1
          else if( cpts2(t)+1 .le. ceiling(ilen(i)+(up-2)*ran(i))) then
            up=max(up-1,1)
          end if
          if(low .ge. up ) then
10022 do 10023 l=0,(low-up)
            help(j+l+up)=0
10023 continue
          end if
          j=j+nint(i)
10024 continue
        end if
        !if(t .ge. 1) return
10028 continue
      end if
      j2=j2-nint(i2)
10029 continue
      return
end subroutine narrowcptsF


recursive subroutine obsnarrowF(n, st, end, le, ri, prev, cu, prevf,res1,res2)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(kind= c_int), intent(in)          :: n, st, end,le,ri
      real(kind= c_double), intent(in)          :: cu(n+1)
      integer(kind= c_int)                      :: g,h,ipm,pl,pr,prev
      real(kind=c_double)                       :: sl,sr,sf,gp,vall,valr,prevf
      integer(kind=c_int), intent(out)          :: res1
      real(kind=c_double), intent(out)          :: res2
      sf=cu(end+1) - cu(st)
      gp=0_c_double
      res1=0
      res2=0_c_double
10030 if(ri-le .gt. 4) goto 10040
10031 do 10032 h=le,min(le+4,end-1)
      sl=cu(h+1)-cu(st)
      sr=cu(end+1)-cu(h+1)
      if( (sl**2 / (h+1-st) + sr**2  / (end-h)) .lt. gp) goto 10032
      gp=sl**2 / (h+1-st) + sr**2  / (end-h)
      ipm=h
10032 continue
      res1=ipm
      res2=sqrt(gp  - sf**2  / (end-st+1))
      goto 10050
10040 continue
      pl=prev-le+1
      pr=ri-prev
      if(pl .gt. pr) then
        g=le+floor(pl/2.0_c_double)
        vall=sqrt( (cu(g+1)-cu(st))**2 / (g+1-st) + (cu(end+1)-cu(g+1))**2 / (end-g) - sf**2 / (end-st+1))
        valr=prevf
      else
        g=prev+floor(pr/2.0_c_double)
        valr=sqrt( (cu(g+1)-cu(st))**2 / (g+1-st) + (cu(end+1)-cu(g+1))**2 / (end-g) - sf**2 / (end-st+1))
        vall=prevf
      endif
      if(vall .le. valr) then
        prevf=valr
        if(pl .gt. pr) then
          call obsnarrowF(n, st, end, le+floor(pl/2.0_c_double), ri, prev, cu, prevf,res1,res2)
        else
          call obsnarrowF(n, st, end, prev, ri, prev+floor(pr/2.0_c_double), cu, prevf,res1,res2)
        endif
      else
        prevf=vall
        if(pl .gt. pr) then
          call obsnarrowF(n, st, end, le, prev, le+floor(pl/2.0_c_double), cu, prevf,res1,res2)
        else
          call obsnarrowF(n, st, end, le, prev+floor(pr/2.0_c_double), prev, cu, prevf,res1,res2)
        endif
      endif
10050 continue
      return
end subroutine obsnarrowF

subroutine sortnarrowF(n, ra, ind)
  use, intrinsic :: ISO_C_BINDING
  implicit none
  !-input/output variables
  integer(kind=c_int), intent(in)   :: n
  integer(kind=c_int), intent(out)  :: ind(n)
  real(kind=c_double)               :: ra(n)
  !-local variables
  real(kind=c_double), parameter :: eps=1e-6
  integer(kind=c_int)               :: i, ir, j, l, iind
  real(kind=c_double)               :: rra
  ! initialize index array
  IF (ind (1) .eq.0) then
     DO i = 1, n
        ind (i) = i
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1
  ir = n
  sorting: do
    ! still in hiring phase
    IF ( l .gt. 1 ) then
       l    = l - 1
       rra  = ra (l)
       iind = ind (l)
       ! in retirement-promotion phase.
    ELSE
       ! clear a space at the end of the array
       rra  = ra (ir)
       !
       iind = ind (ir)
       ! retire the top of the heap into it
       ra (ir) = ra (1)
       !
       ind (ir) = ind (1)
       ! decrease the size of the corporation
       ir = ir - 1
       ! done with the last promotion
       IF ( ir .eq. 1 ) then
          ! the least competent worker at all !
          ra (1)  = rra
          !
          ind (1) = iind
          exit sorting
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l
    ! set up to place rra in its proper level
    j = l + l
    !
    DO while ( j .le. ir )
       IF ( j .lt. ir ) then
          ! compare to better underling
          IF ( ra (j) .gt.  ra (j + 1) )  then
             j = j + 1
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( rra .gt. ra (j) )  then
          ra (i) = ra (j)
          ind (i) = ind (j)
          i = j
          j = j + j
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1
       ENDIF
    ENDDO
    ra (i) = rra
    ind (i) = iind
  END DO sorting
end subroutine sortnarrowF
