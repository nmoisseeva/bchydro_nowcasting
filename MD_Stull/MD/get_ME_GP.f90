      SUBROUTINE get_ME_GP(fstZfile,ni,nj,topo,lat,lon, &
                               latrot,lonrot,AA,xyzgr_rot_8,gdgem)
!*********************************************************************
! author:	Xingxiu Deng
! date:		January 2007
! Object:	to read input topo files, and returns topo,
!               latitude, longitude and rotation matrix for the grid
!               and calculate the coordinates for each grid point
!*********************************************************************
! input:        fstZfile: input FST path and file name
!               ni: No. of grid points along x-dir
!               nj: No. of grid points along y-dir
! output:       topo: filtered topograph of the grid
!               lat: latitude in geographical polar coord
!               lon: longitude in geographical polar coord
!               latrot: latitude in rotated coord
!               lonrot: longitude in rotated coord
!               xyzgr_rot_8: rotated Cartesian coords(3) of grid points
!               gdgem: grid id
! note:
!               may remove lat and lon later as they are not needed from
!               arguments of this subroutine and in MD_prep
!               correspondingly
!*********************************************************************
      implicit none
!*** input
      character(len=100) fstZfile
      integer ni,nj
!*** output
      real topo(ni,nj)
      real lat(ni,nj),lon(ni,nj)  ! geographical latlon
      real latrot(nj),lonrot(ni)  ! rotated latlon
      real*8 AA(3,3)              ! rotation matrix
      real*8 xyzgr_rot_8(ni,nj,3) ! rotated Cartesian coords(3) of GP
      integer gdgem               ! grid id

!*** misc
      integer ix,jy
      real rg12
      logical :: debug=.true.
      integer :: iun=10

!***
      print*, '---- PROCESSING FST FILE --- ', trim(fstZfile)

!*** to get latlon and rotation matrix
      call sub_rotation(gdgem,lat,lon,lonrot,latrot,AA,ni,nj,fstZfile)
      IF (debug) THEN
        print *, ' lonrot(1),latrot(1) = ',lonrot(1),latrot(1)
        print *, ' lonrot(ni),latrot(nj) = ',lonrot(ni),latrot(nj)
        print *, ' lat(1,1),lon(1,1) = ',lat(1,1),lon(1,1)
        print *, ' lat(ni,nj),lon(ni,nj) = ',lat(ni,nj),lon(ni,nj)
      ENDIF

!*** to read ME
      call read_ME_FST_G(fstZfile,iun,ni,nj,topo)

! to calculate each grid point's Cartesian coordinates (xyzgr_rot_8)
! with its origin at the center of the Earth (for GEM grid)
      DO jy=1,nj
      DO ix=1,ni
        call llacar(xyzgr_rot_8(ix,jy,:),lonrot(ix),latrot(jy),1,1)
      ENDDO
      ENDDO

      IF (debug) THEN
        call cal_dist_sphere(xyzgr_rot_8(ni/2,nj/2,:),  &
                             xyzgr_rot_8(ni/2+1,nj/2,:),rg12)
        print*, 'dx=rg12= ', rg12
        call cal_dist_sphere(xyzgr_rot_8(ni/2,nj/2,:),  &
                             xyzgr_rot_8(ni/2,nj/2+1,:),rg12)
        print*, 'dy=rg12=', rg12
      ENDIF

      DO jy=1,nj
      DO ix=1,ni
         IF (topo(ix,jy) .LT. 0.) THEN
           topo(ix,jy)=0.
         ENDIF
      ENDDO
      ENDDO

      IF (debug) THEN
       DO ix=1,3
         print*, (AA(ix,jy),jy=1,3)
       ENDDO
      ENDIF

      RETURN
      END SUBROUTINE get_ME_GP
