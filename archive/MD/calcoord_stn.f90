      SUBROUTINE calcoord_stn(xyz_rot_8,latrot,lonrot,lat,lon,AA)
!*********************************************************************
! author:	Xingxiu Deng
! date:		February 2007
! Object:	To calculate lat and lon in a rotated coordinate defined
!               by rotation matrix given a lat and lon in geographic
!               coordinate system
!*********************************************************************
! input:        lat: latitude of a station
!               lon: longitude of a station
!               AA: rotation matrix
! output:       latrot: latitude of a station in a rotated system
!               lonrot: longitude of a station in a rotated system
!
! intermediate: xyz_8: geographic Cartesian coords(3) of the station
!               xyz_rot_8: rotated Cartesian coords(3) of the station
!               
!*********************************************************************
      implicit none
!*** output
      real latrot, lonrot

!*** input
      real lat, lon
      real*8 AA(3,3)

!*** intermediate
      real*8 xyz_8(3)
      real*8 xyz_rot_8(3)

!*** misc
      logical :: debug=.true.

      print*, '---- Calculating coords for stations -------------'

!--------------------------------------------------------------
      call llacar(xyz_8,lon,lat,1,1)
      call gc2rc(xyz_rot_8,xyz_8,3,AA)
      call rc2rp(latrot,lonrot,xyz_rot_8,3)
      
      if (debug) then
        print*, 'latrot= ', latrot, ' lonrot= ', lonrot
      endif

      RETURN
      end SUBROUTINE calcoord_stn

