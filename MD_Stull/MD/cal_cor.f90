      subroutine cal_cor(Sdaughter,rg,corr,                         &
                     nobsng,Sdaughter_3D,rg_3D,                     &
                     nx,ny,latrot,lonrot,xyzgr_rot_8,topo,Rxy,      &
                     latstnrot,lonstnrot,xyzstnrot_8,topostn,       &
                     a,b,zref1,zref2,avehpbl)
!*********************************************************************
! author:	Xingxiu Deng
! date:		April 2007
! Object:	to calculate sharing factors, circuitous travel distance,
!               and correlation from each of the four coarse source points
!               (around the obs station) to the obs station
!*********************************************************************
! input:        
!       for surface stations:
!               latstnrot: rotated latitude of stations
!               lonstnrot: rotated longitude of stations
!               xyzstnrot_8: rotated Cartesian coords(3) of station
!               topostn: elevation of stations (m)
!       for work grid of MD: 
!               nx,ny: dimension of the fine grid
!               latrot: rotated lat
!               lonrot: rotated lon
!               xyzgr_rot_8: rotated Cartesian coords(3)
!               topo: topography of the fine grid (m)
!
!               Rxy: Radius of influence (km)
!       for MD free parameters:
!               a, b: parameters controling analysis decorrelation rate
!               zref1: terrain-following BL-depth parameter (m)
!               zref2: level-top BL-depth parameter (m) 
!               avehpbl: approximated average BL depth (m)
!       for output from MD_sph:
!               Sdaughter_3D: sharing factors at working GPs from all source points
!               rg_3D: circuitous travel distance from source point to the GP
!               nobsng: No. of source points
! output:       Sdaughter: sharing factor from source point to one obs stn
!               rg: circuitous travel distance from source point to one obs stn
!               corr: correlation from source point to one obs stn
!*********************************************************************
      implicit none
!*** input 
!* - for surface stations
      real latstnrot,lonstnrot
      real topostn
      real*8 xyzstnrot_8(3) ! rotated Cartesian coords(3) of station

!* - for working grid of MD)
      integer nx,ny
      real latrot(ny),lonrot(nx)
      real*8 xyzgr_rot_8(nx,ny,3)
      real topo(nx,ny)

      real Rxy

!* - for output from MD_sph
      integer nobsng
      real :: Sdaughter_3D(nx,ny,nobsng)
      real :: rg_3D(nx,ny,nobsng)
!* free parameters
      real :: a, b, zref1, zref2, avehpbl 

!*** output 
      real :: Sdaughter(nobsng)
      real :: rg(nobsng)
      real :: corr(nobsng)

!***   intermediate
      integer isw_f, jsw_f,ireturnf
      integer iif_index, jjf_index
      integer :: iif_obstn       ! i-index of the obstn's nearest GP within fine grid
      integer :: jjf_obstn       ! j-index of the obstn's nearest GP within fine grid
      real :: Rnearest_obstn     ! distance between the GP at (iif_obstn, jjf_obstn) and obstn
      real :: Znearest_obstn     ! height diff between the GP at (iif_obstn, jjf_obstn) and obstn
      real R_gaussian            ! same as Rxy but in meters

!*** misc
      integer :: isp
      logical :: dbg=.true.

      print*, '---- EXECUTE cal_cor -------------'

      R_gaussian = Rxy * 1.0E03

!*   to find the index of the grid point to the south-west of the
!    station within the fine grid
        CALL find_neighb(isw_f,jsw_f,latstnrot,lonstnrot,      &
                         latrot,lonrot,nx,ny,ireturnf)

!*** Find nearest GP to the station if it is not collcated with a GP
!    based on minimum elevation difference, and get the distance from 
!    the station to the nearest GP and the elevation difference between
!    them
        IF (lonstnrot .NE. lonrot(isw_f) .AND.          &
            latstnrot .NE. latrot(jsw_f) ) THEN
           CALL find_index_HA_sph(isw_f,jsw_f,xyzstnrot_8,topostn, &
                             nx,ny,xyzgr_rot_8,topo,    &
                             iif_index,jjf_index,Rnearest_obstn)
           iif_obstn=iif_index
           jjf_obstn=jjf_index
           Znearest_obstn=topostn-topo(iif_index,jjf_index)
        ELSE
           iif_obstn=isw_f
           jjf_obstn=jsw_f
           Rnearest_obstn=0.0
           Znearest_obstn=0.0
        ENDIF
        
!*** calculate sharing factors from each of the four coarse source points 
!    (around the obs station) to the obs station
        DO isp=1,nobsng
          Sdaughter(isp)=Sdaughter_3D(iif_obstn,jjf_obstn,isp)*        &
             exp(-0.5*Znearest_obstn*Znearest_obstn/(avehpbl*avehpbl))
          rg(isp)=rg_3D(iif_obstn,jjf_obstn,isp)+Rnearest_obstn
        ENDDO

!*** calculate correlation from each of the four coarse source points 
!    (around the obs station) to the obs station
        DO isp=1,nobsng
          corr(isp) = Sdaughter(isp)*     &
             exp(-0.5*rg(isp)*rg(isp)/(R_gaussian*R_gaussian))
        ENDDO

      RETURN
      END
