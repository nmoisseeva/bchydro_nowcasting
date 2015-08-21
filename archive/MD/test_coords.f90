      SUBROUTINE test_coords(indir,nxc,nyc,xyzgrc_rot_8,AAc,nxf,nyf, &
                             xyzgrf_rot_8,AAf,AA)
!*********************************************************************
! author:	Xingxiu Deng
! date:		January 2007
! Object:	to read in gfactor from the output file of design_grid
!               and test if the fine grid is overlapped with the coarse grid
!*********************************************************************
! input:        indir: directory where gc2f_ratio is stored
!               nxc, nyc: dimension size for the coarse grid
!               xyzgrc_rot_8: rotated Cartesian coords for the coarse grid
!               nxf, nyf: dimension size for the fine grid
!               xyzgrf_rot_8: rotated Cartesian coords for the fine grid
! derived input:
!               gfactor: grid spacing ratio of the coarse to fine grid
!                        it is set in MD.cfg and is output from
!                        design_grid.ksh
!*********************************************************************
      implicit none
!*** input
      character(len=70) indir
      integer nxc,nyc
      real*8 xyzgrc_rot_8(nxc,nyc,3) 
      real*8 AAc(3,3)
      integer nxf,nyf
      real*8 xyzgrf_rot_8(nxf,nyf,3) 
      real*8 AAf(3,3)

      real*8 AA(3,3)

!*** derived from input: grid spacing ratio of the coarse to fine grid
      integer gfactor

!*** misc
      integer k,iic,jjc,iif,jjf
      real dif1, dif2, dif3
      real :: epl=1.0E-4
      integer i,j

      open(21,file=trim(indir)//'/gc2f_ratio',form='formatted',   &
              status='old')
      read(21,'(i)') gfactor 
      close(21)
      print*, 'gfactor= ', gfactor

      iic=nxc
      jjc=nyc
      iif=(iic-1)*gfactor+1
      jjf=(jjc-1)*gfactor+1
      dif1=xyzgrc_rot_8(iic,jjc,1) - xyzgrf_rot_8(iif,jjf,1)
      dif2=xyzgrc_rot_8(iic,jjc,2) - xyzgrf_rot_8(iif,jjf,2)
      dif3=xyzgrc_rot_8(iic,jjc,3) - xyzgrf_rot_8(iif,jjf,3)
      IF (dif1 .GT. epl .OR. dif2 .GT. epl .OR. dif3 .GT. epl) THEN
        print*, 'the fine grid is not overlapped with the coarse grid'
        print*, 'iic,jjc= ', iic,jjc,(xyzgrc_rot_8(iic,jjc,k),k=1,3)
        print*, 'iif,jjf= ', iif,jjf,(xyzgrf_rot_8(iif,jjf,k),k=1,3)
        print*, 'check the diff again to see if it is OK with small epl'
        stop
      ENDIF
      
      DO j=1,3
      DO i=1,3
        IF ( AAc(i,j) .NE. AAf(i,j) ) THEN
         print*,'The fine grid should have same rotation to the coarse'
         print*, 'check design_grid'
         stop
        ELSE
         AA(i,j)=AAf(i,j)
        ENDIF
      ENDDO
      ENDDO

      RETURN
      end
