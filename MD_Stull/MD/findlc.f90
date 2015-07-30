SUBROUTINE findlc(nx,ny,xgr,ygr,xpt,ypt,ipt,jpt,ireturn)
!*********************************************************************
! author:	Xingxiu Deng
! date:		January 2003
! Object:	to find index of the grid point to the south-west of
!               a given station within a grid
!*********************************************************************
! input:        nx, ny: dimension size of the grid
!               xgr: 2D x coordinate of grid points
!               ygr: 2D y coordinate of grid points
!               xpt: x coord of a station
!               ypt: y coord of a station
! output:       ipt: i-index in the grid to the south-west of the station
!               jpt: j-index in the grid to the south-west of the station
!               ireturn: < 0, the station is outside of the domain
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
      REAL :: xpt
      REAL :: ypt
!*** output
      INTEGER :: ipt
      INTEGER :: jpt
      INTEGER :: ireturn
!*** misc 
      INTEGER :: i,j
      INTEGER :: ipt_1D(ny), jpt_1D(nx)

      ireturn = -1
      ipt_1D = 0
      jpt_1D = 0

      DO j=1, ny-1
      DO i=1, nx-1
         IF (xpt .GE. xgr(i,j) .AND. xpt .LE.xgr(nx-1,j)) THEN
            ipt_1D(j) =i
         ENDIF
      ENDDO
      ENDDO

      DO j=1, ny-1
       IF(ypt .GE. ygr(ipt_1D(j),j) .AND. ipt_1D(j) .GT. 0  &
          .AND. ypt .LE. ygr(ipt_1D(j),ny-1)) THEN
         ipt = ipt_1D(j)
       ENDIF
      ENDDO
      DO i=1, nx-1
      DO j=1, ny-1
        IF (ypt .GE. ygr(i,j).AND. ypt .LE.ygr(i,ny-1) ) THEN
          jpt_1D(i) =j
        ENDIF
      ENDDO
      ENDDO
      DO i=1, nx-1
        IF(xpt .GE. xgr(i,jpt_1D(i)) .AND. jpt_1D(i).GT. 0  &
           .AND. xpt .LE. xgr(nx-1,jpt_1D(i)) ) THEN
         jpt = jpt_1D(i)
         ireturn  = 0
        ENDIF
      ENDDO

      RETURN
    END SUBROUTINE findlc
