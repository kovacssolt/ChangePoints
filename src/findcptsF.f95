subroutine findcptsF(x,n,dec,dep,ilen,nint,nsum,par,stats,pen,thr,t,t2,tot,bou,bou2,cpts,cpts2) bind(c,name="findcptsF_")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      real(kind=c_double), parameter :: eps=1e-4
      integer(kind= c_int), intent(in)          :: n,dep,nsum,par,pen,stats
      integer(kind= c_int), intent(in)          :: nint(dep-1)
      real(kind=c_double),  intent(in)          :: x(n),dec,ilen(dep-1),thr
      integer(kind= c_int)                      :: i,j,k,l,m,ipm,prev,res1,low,up,j2,t,t2
      real(kind=c_double)                       :: cu(n+1),ran(dep-1),temp(n-1),sl,sr,sf,prevf,res2,n2
      integer(kind= c_int),  dimension (:), allocatable :: ord,s(:,:)
      real(kind= c_double), intent(out)         :: bou(nsum),cpts(n,3),tot
      integer(kind= c_int), intent(out)         :: bou2(nsum,3),cpts2(n,3)
10000 j=2
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
        call obsF(n, bou2(1,1), bou2(1,2), bou2(1,1), bou2(1,2), prev, cu, prevf,res1,res2)
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
10007 do 10014 i=1,dep-1
      ran(i) = (n-ilen(i))/(nint(i)-1)
      if(par .lt. ilen(i)+1) then
10008 do 10009 k=1,nint(i)
        bou2(j,1)=floor(1+(k-1)*ran(i))
        bou2(j,2)=ceiling(ilen(i)+(k-1)*ran(i))
        sf=cu(bou2(j,2)+1) - cu(bou2(j,1))
        prev=bou2(j,1)+floor((bou2(j,2)-bou2(j,1))/3.0_c_double)
        sl=cu(prev+1)-cu(bou2(j,1))
        sr=cu(bou2(j,2)+1)-cu(prev+1)
        prevf=sqrt( sl**2 / (prev+1-bou2(j,1)) + sr**2 / (bou2(j,2)-prev) - sf**2 / (bou2(j,2)-bou2(j,1)+1) )
        call obsF(n, bou2(j,1), bou2(j,2), bou2(j,1), bou2(j,2), prev, cu, prevf,res1,res2)
        bou2(j,3)=res1
        bou(j)=res2
        j=j+1
10009 continue
      else
10010 do 10013 k=1,nint(i)
        bou2(j,1)=floor(1+(k-1)*ran(i))
        bou2(j,2)=ceiling(ilen(i)+(k-1)*ran(i))
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
        j=j+1
10013 continue
      endif
10014 continue
      if(stats .eq. 1) return
      cpts=0.0_c_double
      cpts2=0
      t=1
10016 call findcpF(n, 1, n, t, dep, ilen, nint, nsum, ran, thr, bou, bou2, cpts, cpts2)
      t=t-1
      tot=sum(x**2)-cu(n+1)**2/n
      if(t .eq. 0) return
      allocate(ord(t))
10017 do 10018 i=1,t
      ord(i)=i
10018 continue
      call sortF(t,cpts(1:t,1),ord)
      sl=cu(cpts2(1,1)+1)
      sr=cu(n+1)-cu(cpts2(1,1)+1)
      cpts2(1:t,1:3)=cpts2(ord,1:3)
      !cpts(1,2)=sum(x**2) -sl**2/cpts2(1,1) -sr**2/(n-cpts2(1,1))
      cpts(1,2)=tot-sl**2/cpts2(1,1) -sr**2/(n-cpts2(1,1)) +cu(n+1)**2/n
10019 do 10020 i=2,t
      sf=cu(cpts2(i,3)+1)-cu(cpts2(i,2))
      sl=cu(cpts2(i,1)+1)-cu(cpts2(i,2))
      sr=cu(cpts2(i,3)+1)-cu(cpts2(i,1)+1)
      cpts(i,2)=cpts(i-1,2) - sl**2/(cpts2(i,1)+1-cpts2(i,2)) - sr**2/(cpts2(i,3)-cpts2(i,1)) + sf**2/(cpts2(i,3)-cpts2(i,2)+1)
10020 continue
      n2=real(n)
      if(pen .eq. 1) then
10021 do 10022 i=1,t
        cpts(i,3)=cpts(i,2)+i*log(n2)
10022 continue
      else if(pen .eq. 2) then
      prevf=0.0_c_double
10023 do 10024 i=1,t
        prevf=prevf-log((cpts2(i,3)-cpts2(i,2)+1)/n2)+log((cpts2(i,1)-cpts2(i,2)+1)/n2)+log((cpts2(i,3)-cpts2(i,1))/n2)
        cpts(i,3)=cpts(i,2)+i*log(n2)*3/2+prevf/2
10024 continue
      end if
      t2=minloc(cpts(1:t,3),1)
      if(cpts(t2,3) .ge. tot) t2=0
      return
end subroutine findcptsF


recursive subroutine obsF(n, st, end, le, ri, prev, cu, prevf,res1,res2)
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
          call obsF(n, st, end, le+floor(pl/2.0_c_double), ri, prev, cu, prevf,res1,res2)
        else
          call obsF(n, st, end, prev, ri, prev+floor(pr/2.0_c_double), cu, prevf,res1,res2)
        endif
      else
        prevf=vall
        if(pl .gt. pr) then
          call obsF(n, st, end, le, prev, le+floor(pl/2.0_c_double), cu, prevf,res1,res2)
        else
          call obsF(n, st, end, le, prev+floor(pr/2.0_c_double), prev, cu, prevf,res1,res2)
        endif
      endif
10050 continue
      return
end subroutine obsF


recursive subroutine findcpF(n, st, end, t, dep, ilen, nint, nsum, ran, thr, bou, bou2, cpts, cpts2)
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(kind= c_int), intent(in)          :: n, st, end, dep, nsum
      integer(kind= c_int), intent(in)          :: nint(dep-1),bou2(nsum,3)
      real(kind= c_double), intent(in)          :: ilen(dep-1), ran(dep-1), bou(nsum), thr
      integer(kind= c_int)                      :: i,j,l,lay(dep),low,up,help,cp
      real(kind= c_double)                      :: lay2(dep)
      integer(kind=c_int), intent(inout)        :: cpts2(n,3),t
      real(kind=c_double), intent(inout)        :: cpts(n,2)
      lay2=0.0_c_double
      j=nsum
      if( (st .eq. 1 ) .and. (end .eq. n)) then
          lay(1)=1
          lay2(1)=bou(1)
      end if
      i=dep-1
10061 subsets: do while ( (i .ge. 1) .AND. (ilen(i) .lt. (end-st+1)) )
      j=j-nint(i)
      low=ceiling((st-1)/ran(i)+1)
      if( st .gt. bou2(j+low,1)) then
      !if( st .gt. floor(1+(low-1)*ran(i))) then
         !low=min(low+1,nint(i)+1)
         low=low+1
      else
        if( st .le. floor(1+(low-2)*ran(i))) then
          low=low-1
        end if
      end if
      up=floor((end-ilen(i))/ran(i)+1)
      if( end .lt. bou2(j+up,2)) then
      !if( end .lt. ceiling(ilen(i)+(up-1)*ran(i))) then
         !up=max(up-1,0)
         up=up-1
      else
        if( end .ge. ceiling(ilen(i)+up*ran(i))) then
        up=up+1
        end if
      end if
      i=i+1
      if(up .ge. low ) then
        lay(i)=maxloc(bou((j+low):(j+up)),1)+j+low-1
        lay2(i)=bou(lay(i))
!10062 do 10064 l=0,(up-low)
!         if(bou(j+l+low) .gt. lay2(i)) then
!           lay(i)=j+l+low
!            lay2(i)=bou(lay(i))
!          end if
!10064 continue
      end if
      i=i-2
      end do subsets
10066 continue
      !help=maxloc(lay2,1)
      help=maxloc(lay2(i+1:dep),1)+i
      if(lay2(help) .le. thr) return
      !if(lay2(help) .eq. 0.0_c_double) return
      cp=bou2(lay(help),3)
      cpts(t,1)=lay2(help)
      cpts2(t,1)=cp
      cpts2(t,2)=st
      cpts2(t,3)=end
      !t einbauen laufvariable
      t=t+1
      call findcpF(n, st, cp, t, dep, ilen, nint, nsum, ran, thr, bou, bou2, cpts, cpts2)
      call findcpF(n, cp+1, end, t, dep, ilen, nint, nsum, ran, thr, bou, bou2, cpts, cpts2)
      return
end subroutine findcpF


subroutine sortF(n, ra, ind)
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
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1
       ENDIF
    ENDDO
    ra (i) = rra
    ind (i) = iind
  END DO sorting
end subroutine sortF
