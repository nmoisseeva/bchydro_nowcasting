       subroutine read_ME_FST_G(fstZfile,iun,nx,ny,ME)
!*********************************************************************
! author:	Xingxiu Deng
! date:		February 2006
! Object:	To open fst Z grid file and read ">>" , "^^", and ME
!*********************************************************************
!Input variables:
!      fstZfile: input FST path and file name
!      iun: unit number
!      nx,ny: dimension size
!---------------------------------------------------------------------
!Output variables:
!      ME: topography
!*********************************************************************
      IMPLICIT none
      include 'fst_variables.cdk'

      integer :: nx,ny
      integer :: ni,nj

      character (len=100) fstZfile
      integer :: iun

      real ME(nx,ny)

      print*, '------- READ ME ---------'

!*     open FST Z grid file 
      print *,"-- Opening FST Z grid File --",trim(fstZfile)
      ier=fnom(iun,trim(fstZfile),'RND',0)

      if (ier.lt.0) then
        print*,'Fatal error while opening the file: ', trim(fstZfile)
        stop
      endif

      nrecs = fstouv (iun, 'RND')
      if (nrecs.lt.1) then
         print*,'Error, no record found in file: ', trim(fstZfile)
         ier = fstfrm(iun)
         ier = fclos(iun)
         stop
      endif

       keyX = fstinf(iun,nc,nr,nl,-1,' ',-1,-1,-1,' ','>>')
       ni=nc

       if (ni .ne. nx) then
         print*,'Dimension size in file ', trim(fstZfile),  &
                ' does not match input nx'
         ier = fstfrm(iun)
         ier = fclos(iun)
         stop
       endif
  
       keyY = fstinf(iun,nc,nr,nl,-1,' ',-1,-1,-1,' ','^^')
       nj=nr
       if (nj .ne. ny) then
         print*,'Dimension size in file ', trim(fstZfile),  &
                ' does not match input ny'
         ier = fstfrm(iun)
         ier = fclos(iun)
         stop
       endif

       ier=fstlir(ME,iun,nc,nr,nl,-1,' ',-1,-1,-1,' ','ME')
       print*, 'ME(1,1)= ', ME(1,1)

! close fst file
 999   CONTINUE
       ier = fstfrm(iun)
       ier = fclos(iun)
       RETURN
       END
