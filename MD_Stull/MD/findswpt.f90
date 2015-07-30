      SUBROUTINE findswpt(nx,ny,xgr,ygr,xstn,ystn,isw,jsw,ireturn)
!*********************************************************************
! author:	Xingxiu Deng
! date:		February 2007
! Object:	 Search in x and y to find isw and jsw of xstn, ystn
!*********************************************************************
! input:        xgr, ygr: 2D x and y coords of the grid  in comp space
!               nx, ny: dimension size of the grid
!               xtsn, ystn: x and y coords of one station
! output:       isw: i-index in the grid to the south-west of the station
!               jsw: j-index in the grid to the south-west of the station
!               ireturn:  =0 if the station is inside the domain
!                         <0 if the station is outside of the domain
!
! Note: x* and y* should be monotonically increasing as i and j
!               indices increase
!               
!*********************************************************************
      IMPLICIT NONE

!*** input
      INTEGER :: nx,ny
      REAL :: xgr(nx,ny)
      REAL :: ygr(nx,ny)
      REAL :: xstn
      REAL :: ystn

!*** output
      INTEGER :: isw
      INTEGER :: jsw
      INTEGER :: ireturn

!*** misc
      INTEGER :: i,j
      INTEGER :: isw_1D(ny), jsw_1D(nx)
!
      ireturn = -1
      isw_1D = 0
      jsw_1D = 0

      DO j=1, ny-1
      DO i=1, nx-1
        IF (xstn .GE. xgr(i,j) .AND. xstn .LE.xgr(nx-1,j)) THEN
          isw_1D(j) =i
        ENDIF
      ENDDO
      ENDDO

      DO j=1, ny-1
        IF(ystn .GE. ygr(isw_1D(j),j) .AND. isw_1D(j) .GT. 0  &
           .AND. ystn .LE. ygr(isw_1D(j),ny-1)) THEN
          isw = isw_1D(j)
        ENDIF
      ENDDO
      DO i=1, nx-1
      DO j=1, ny-1
        IF (ystn .GE. ygr(i,j).AND. ystn .LE.ygr(i,ny-1) ) THEN
          jsw_1D(i) =j
        ENDIF
      ENDDO
      ENDDO
      DO i=1, nx-1
        IF(xstn .GE. xgr(i,jsw_1D(i)) .AND. jsw_1D(i).GT. 0  &
           .AND. xstn .LE. xgr(nx-1,jsw_1D(i)) ) THEN
         jsw = jsw_1D(i)
         ireturn  = 0
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE findswpt
