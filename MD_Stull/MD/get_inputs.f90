      SUBROUTINE get_inputs(outfstpf,topopf_coarse,topopf_fine,sfcstapf, &
                            maxsta,Rxy,indir)
!*********************************************************************
! author:	Xingxiu Deng
! date:		January 2007
! Object:	to read input topo files for the coarse and fine grids,
!               to calculate their rotated latlon and cartesian coordinates,
!               to test whether or not the fine grid is overlapped with the
!               coarse grid
!               and to read the list file for sfc stations
!*********************************************************************
! input:        topopf_coarse: path/file name for coarse grid topo file
!               topopf_fine: path/file name for fine grid topo file
!               sfcstapf: path/file name containing the list of all sfc
!                           stations
!               Rxy: Radius of influence (km)
!               indir: directory where holds inputs
! output:       outfstpf: path/file name of output sharing factor file in
!                         RPN standard format
! Derived from input:
!    coarse grid:
!               nxc,nyc: dimension size for the coarse grid
!               topoc: filtered topography for the coarse grid
!               latc, lonc: latlon for the coarse grid
!               latrotc, lonrotc: rotated latlon for the coarse grid
!               xyzgrc_rot_8: rotated Cartesian coords for the coarse grid
!               AAc: rotation matrix for the coarse grid
!    fine grid:
!               nxf,nyf: dimension size for the fine grid
!               topof: filtered topography for the fine grid
!               latf, lonf: latlon for the fine grid
!               latrotf, lonrotf: rotated latlon for the fine grid
!               xyzgrf_rot_8: rotated Cartesian coords for the fine grid
!               AAf: rotation matrix for the fine grid
!    surface station:
!        with maximum dimension size of maxstn
!               latstn_m: latitude of stations
!               lonstn_m: longitude of stations
!               elevstn_m: elevation of stations
!               stname_m: name or character id of stations
!               nobs: No. of surface stations
!        with true dimension size of nobs
!               latstn: latitude of stations
!               lonstn: longitude of stations
!               elevstn: elevation of stations
!               stname: name or character id of stations
!               
!*********************************************************************
      implicit none
!*** output
      character(len=100) outfstpf

!*** input
      character(len=100) topopf_coarse, topopf_fine
      character(len=100) sfcstapf
      integer maxsta
      real Rxy
      character(len=70) indir

!*** derived from input
!** for coarse grid
      integer nxc,nyc
      real,allocatable :: topoc(:,:)
      real,allocatable :: latc(:,:),lonc(:,:) ! geographical polar coord
      real*8,allocatable :: xyzgrc_8(:,:,:)   ! geographical Cartesian coords(3) of GP
      real,allocatable :: latrotc(:),lonrotc(:) ! rotated polar coord
      real*8,allocatable :: xyzgrc_rot_8(:,:,:) ! rotated Cartesian coords(3) of GP
      real*8 AAc(3,3)                           ! rotation matrix
      integer gdgem_c                           ! grid id 

!** for fine grid
      integer nxf,nyf
      real,allocatable :: topof(:,:)
      real,allocatable :: latf(:,:),lonf(:,:) ! geographical polar coord
      real*8,allocatable :: xyzgrf_8(:,:,:)   ! geographical Cartesian coords(3) of GP
      real,allocatable :: latrotf(:),lonrotf(:) ! rotated polar coord
      real*8,allocatable :: xyzgrf_rot_8(:,:,:) ! rotated Cartesian coords(3) of GP
      real*8 AAf(3,3)                           ! rotation matrix
      integer gdgem_f                           ! grid id 

!** for coarse and fine grid
      real*8 AA(3,3)                            ! rotation matrix shared by the coarse and fine grid

!** for surface stations
      integer nobs
      real, allocatable :: latstn_m(:),lonstn_m(:),elevstn_m(:)
      character (LEN=8), allocatable :: stname_m(:)

      real,allocatable :: latstn(:),lonstn(:),elevstn(:)
      character (LEN=8), allocatable :: stname(:)

!*** misc
      integer iun
      parameter (iun=20)
      integer ista,k

      print*, '---- EXECUTE get_inputs -------------'

!--------------------------------------------------------------
!*** get dimension size for coarse grid
      call getninj_Z(topopf_coarse,iun,nxc,nyc)

!*** allocate memory for coarse grid
      allocate(topoc(nxc,nyc))
      allocate(latc(nxc,nyc))
      allocate(lonc(nxc,nyc))
      allocate(xyzgrc_8(nxc,nyc,3))
      allocate(latrotc(nyc))
      allocate(lonrotc(nxc))
      allocate(xyzgrc_rot_8(nxc,nyc,3))

!*** get dimension size for fine grid
      call getninj_Z(topopf_fine,iun,nxf,nyf)

!*** allocate memory for fine grid
      allocate(topof(nxf,nyf))
      allocate(latf(nxf,nyf))
      allocate(lonf(nxf,nyf))
      allocate(xyzgrf_8(nxf,nyf,3))
      allocate(latrotf(nyf))
      allocate(lonrotf(nxf))
      allocate(xyzgrf_rot_8(nxf,nyf,3))

!--------------------------------------------------------------
!*** read the coarse grid topo file and return ME, latlonrot,
!*** rotation matrix, and rotated Cartesian coords

!**  for the coarse grid
      call get_ME_GP(topopf_coarse,nxc,nyc,topoc, &
                     latc,lonc,latrotc,lonrotc,AAc,xyzgrc_rot_8,gdgem_c)

!**  for the fine grid
      call get_ME_GP(topopf_fine,nxf,nyf,topof, &
                     latf,lonf,latrotf,lonrotf,AAf,xyzgrf_rot_8,gdgem_f)

!*** to test whether or not the fine grid is overlapped with the coarse grid
      call test_coords(indir,nxc,nyc,xyzgrc_rot_8,AAc,nxf,nyf,  &
                       xyzgrf_rot_8,AAf,AA)

!--------------------------------------------------------------
!*** read the list of surface observation stations
      print*, '---- READ SURFACE OBSERVATION STATIONS -------------'
      allocate(stname_m(maxsta))
      allocate(latstn_m(maxsta))
      allocate(lonstn_m(maxsta))
      allocate(elevstn_m(maxsta))

      call read_sfcstation(sfcstapf,maxsta,nobs,stname_m,latstn_m, &
                           lonstn_m,elevstn_m)

      allocate(stname(nobs))
      allocate(latstn(nobs))
      allocate(lonstn(nobs))
      allocate(elevstn(nobs))

      DO ista=1,nobs
         stname(ista)=stname_m(ista)
         latstn(ista)=latstn_m(ista)
         lonstn(ista)=lonstn_m(ista)
         elevstn(ista)=elevstn_m(ista)
      ENDDO

      deallocate(stname_m)
      deallocate(latstn_m)
      deallocate(lonstn_m)
      deallocate(elevstn_m)

!--------------------------------------------------------------
!*** To find neighbouring coarse grid points around each surface
!    station ad treat them as source points to prepare for MD
      call source_pt(outfstpf,nobs,stname,latstn,lonstn,elevstn,  &
                   nxc,nyc,latrotc,lonrotc,xyzgrc_rot_8,topoc,gdgem_c, &
                   nxf,nyf,latrotf,lonrotf,xyzgrf_rot_8,topof,gdgem_f, &
                   AA,Rxy,indir)

      RETURN
      end SUBROUTINE get_inputs
