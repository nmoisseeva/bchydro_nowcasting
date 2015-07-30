      SUBROUTINE source_pt(outfstpf,nobs,stname,latstn,lonstn,elevstn,  &
                   nxc,nyc,latrotc,lonrotc,xyzgrc_rot_8,topoc,gdgem_c,&
                   nxf,nyf,latrotf,lonrotf,xyzgrf_rot_8,topof,gdgem_f,&
                   AA,Rxy,indir)
!*********************************************************************
! author:	Xingxiu Deng
! date:		February 2007
! Object:	To find neighbouring coarse grid points around each 
!               surface station ad call SUB MD_prep to prepare for MD
!*********************************************************************
! input:        
!       for surface stations:
!               nobs: No. of surface stations
!               stname: name or character id of stations
!               latstn: latitude of stations
!               lonstn: longitude of stations
!               elevstn: elevation of stations (m)
!       derived for surface stations
!               latstn_rot: rotated lat of stn within coarse/fine grid
!               lonstn_rot: rotated lon of stn within coarse/fine grid
!               xyzstn_rot_8:  rotated Cartesian coords(3) of station
!               iswc: i index of the GP to the south and west of stn (wrt coarse grid)
!               jswc: j index of the GP to the south and west of stn (wrt coarse grid)
!               ireturnc: < 0, the station is outside the coarse grid 
!                         iswc=-999 and jswc=-999.
!       for coarse grid: see below
!       for fine grid: see below
!
!               Rxy: Radius of influence (km)
!               indir: directory where holds inputs
! output:       outfstpf: path/file name of output sharing factor file in
!                         RPN standard format
!*********************************************************************
      implicit none
!*** output
      character(len=100) outfstpf

!*** input
      real Rxy
      character(len=70) indir

!** for coarse grid
      integer nxc,nyc
      real topoc(nxc,nyc)             ! topography (m)
      real latrotc(nyc),lonrotc(nxc)  ! rotated polar coord
      real*8 xyzgrc_rot_8(nxc,nyc,3)  ! rotated Cartesian coords(3) of GP
      integer gdgem_c                 ! grid id

!** for fine grid
      integer nxf,nyf
      real topof(nxc,nyc)             ! topography (m)
      real latrotf(nyc),lonrotf(nxc)  ! rotated polar coord
      real*8 xyzgrf_rot_8(nxc,nyc,3)  ! rotated Cartesian coords(3) of GP
      integer gdgem_f                 ! grid id

!** for coarse and fine grid
      real*8 AA(3,3)                 ! rotation matrix shared by the coarse and fine grid

!** for surface stations
      integer nobs
      real :: latstn(nobs),lonstn(nobs),elevstn(nobs)
      character (LEN=8) :: stname(nobs)

!** derived for surface stations
      real latstn_rot(nobs),lonstn_rot(nobs)
      real*8 xyzstn_rot_8(nobs,3) ! rotated Cartesian coords(3) of station
      integer iswc(nobs),jswc(nobs)
      integer ireturnc(nobs)

!** functions
      integer gdxyfll,ier

!*** misc
      integer ista

      print*, '---- EXECUTE MD_prep -------------'

!--------------------------------------------------------------
      do ista=1,nobs
         print*, '--- PROCESSING STATION ', ista, ' - ', stname(ista)  

!*** to calculate each station's coordinates within the rotated grid
         call calcoord_stn(xyzstn_rot_8(ista,:),latstn_rot(ista),     &
                          lonstn_rot(ista),latstn(ista),lonstn(ista),AA)

!---------------------------------------------------------------------
!*** gdxyfll returns x and y coord in grid point unit, rather than
!*   rotated lat/lon, we can use it only if change SUB find_neighb and
!*   all rest of SUBROUTINES that use rotated lat/lon, no change for now
!
!         ier=gdxyfll(gdgem_c,lonstn_rot(ista),latstn_rot(ista),   &
!                     latstn(ista),lonstn(ista),1)
!         print*, 'wrt coarse grid ', 'latrot= ', latstn_rot(ista), &
!                 ' lonrot= ', lonstn_rot(ista)
!---------------------------------------------------------------------

!*** to find the index of the grid point to the south-west of the
!    station within the coarse grid
         call find_neighb(iswc(ista),jswc(ista),latstn_rot(ista),  &
                          lonstn_rot(ista),latrotc,lonrotc,nxc,nyc,&
                          ireturnc(ista))
         print*, 'iswc, jswc= ', iswc(ista),jswc(ista), ireturnc(ista)
      enddo

!*** to prepare for calling the main program MD
       call MD_prep(outfstpf,nobs,stname,latstn_rot,lonstn_rot,elevstn, &
                    xyzstn_rot_8,iswc,jswc,ireturnc,nxc,nyc,topoc,nxf,nyf, &
                    topof,latrotf,lonrotf,xyzgrf_rot_8,Rxy,indir)

      RETURN
      end
