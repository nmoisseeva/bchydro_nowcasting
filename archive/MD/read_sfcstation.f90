      SUBROUTINE read_sfcstation(infile,maxsta,nobs,stn,lat,lon,elev)
!*********************************************************************
! author:	Xingxiu Deng
! date:		December 14, 2006
!
! Object:	to read the list of surface stations
!*********************************************************************
! input:        infile: input file for the list of surface stations
!               maxsta: maximum number of stations allowed
!
! output:       nobs: total number of surface stations
!               stn:  station name or id 
!               lat:  station latitude (-90 - 90)
!               lon:  station longitude (-180 - 180)
!               elev: station elevation (meter)
!*********************************************************************
      
      IMPLICIT NONE
      CHARACTER (LEN=100) :: infile
      INTEGER :: maxsta
      REAL :: lat(maxsta),lon(maxsta),elev(maxsta)
      CHARACTER (LEN=8) :: stn(maxsta)
      INTEGER :: nobs
      INTEGER :: istat
      INTEGER :: ierror, k
      LOGICAL :: debug = .true.

      istat = 0

      OPEN(10,file=trim(infile),form='formatted', &
              status='old',IOSTAT=ierror)
      IF (ierror.ne.0) THEN
       print*, 'could not open file: ', infile
       print*, 'system error code: ',ierror
       istat = -1
       RETURN
      ENDIF

      READ(10,800) nobs
 800  FORMAT(1x, i8)
      print*,'There are ', nobs, ' surface stations'

      IF (nobs > maxsta) THEN
       print 890, maxsta,nobs
 890   FORMAT('** error in read_sfcstation: maxsta = ',i8,/, &
              'but there are ',i8,' stations in the input file',/)
       print*,'** increase maxsta and try again'
       istat = -2
       RETURN
      END IF

      DO k=1,nobs
      READ(10,*) stn(k),lat(k),lon(k),elev(k)
! 801  FORMAT(1X,a8,f10.5,1X,f10.4,1X,f10.4) 
      ENDDO

      IF ( debug ) THEN
         DO k=1,nobs
          print*, stn(k),lat(k),lon(k),elev(k)
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE read_sfcstation
