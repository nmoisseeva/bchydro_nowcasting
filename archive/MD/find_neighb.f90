      SUBROUTINE find_neighb(isw,jsw,latrotstn,lonrotstn,latrot,  &
                             lonrot,ni,nj,ireturn)
!*********************************************************************
! author:	Xingxiu Deng
! date:		February 2007
! Object:	to find index of the grid point to the south-west of
!               a given station within a grid
!*********************************************************************
! input:        latrotstn: latitude of a station in a rotated system
!               lonrotstn: longitude of a station in a rotated system
!               latrot: latitude of a rotated grid
!               lonrot: longitude of a rotated grid
!               ni, nj: dimension size of the grid
! output:       isw: i-index in the grid to the south-west of the station
!               jsw: j-index in the grid to the south-west of the station
!               ireturn: < 0, the station is outside of the domain, no
!                        isw or jsw is returned.
!
! intermediate: xs: x coord of a station
!               ys: y coord of a station
!               xgr: 2D x coordinate of grid points
!               ygr: 2D y coordinate of grid points
!
! Note: x* and y* should be monotonically increasing as i and j
!               indices increase
!               
!*********************************************************************
      implicit none
!*** output
      integer isw,jsw,ireturn

!*** input
      integer ni,nj
      real latrotstn,lonrotstn
      real latrot(nj),lonrot(ni)

!*** intermediate
      real xs, ys
      real xgr(ni,nj), ygr(ni,nj)
      REAL :: xgrmin,xgrmax,ygrmin,ygrmax

!*** misc
      logical :: debug=.true.
      integer i,j

      print*, '---- find_neighb  -------------'

!* define x* and y* coords
!* since the input lat/lon is wrt rotated computational spac, so
!* directly use latrot and lonrot for x* and y*
      xs=lonrotstn
      ys=latrotstn
      DO j=1,nj
      DO i=1,ni
        xgr(i,j)=lonrot(i)
        ygr(i,j)=latrot(j)
      ENDDO
      ENDDO

!* find xgrmin,xgrmax,ygrmin,ygrmax
      xgrmin = 1.0E+8
      xgrmax = -1.0E+8
      ygrmin = 1.0E+8
      ygrmax = -1.0E+8

      DO j =1,nj
        IF (xgr(1,j) .LT. xgrmin) xgrmin =  xgr(1,j)
        IF (xgr(ni,j) .GT. xgrmax) xgrmax =  xgr(ni,j)
      ENDDO      
      DO i =1,ni
        IF (ygr(i,1) .LT. ygrmin) ygrmin =  ygr(i,1)
        IF (ygr(i,nj) .GT. ygrmax) ygrmax =  ygr(i,nj)
      ENDDO  

!* find grid point to the south and west of the station 
      IF (xs .LT. xgrmin .OR. xs .GT. xgrmax .OR.    &
          ys .LT. ygrmin .OR. ys .GT. ygrmax  ) THEN
          ireturn = -1
          print*,'xgrmin,xgrmax,ygrmin,ygrmax',xgrmin,xgrmax,ygrmin,ygrmax
          print*,'xs,ys',xs,ys
          isw=-999
          jsw=-999
      ELSE
          CALL findswpt(ni,nj,xgr,ygr,xs,ys,isw,jsw,ireturn)
      ENDIF

      if (debug) then
        print*, 'isw= ', isw, ' jsw= ', jsw,' ireturn= ', ireturn
      endif

      RETURN
      end SUBROUTINE find_neighb
