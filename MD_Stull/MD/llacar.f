***s/r llacar - transformation from a set of points (F_lat,F_lon) in 
*             the spherical coordinate system to cartesian space
*
      subroutine llacar( F_xyz_8,F_lon,F_lat,ni,nj)
      implicit none
      integer ni,nj
      real*8 F_xyz_8(3,ni*nj)
      real F_lon(ni),F_lat(nj)
*
*author 
*     Michel Roch - April 90
*
*revision
* v2_00 - Lee V.            - initial MPI version (from llacar v1_03)
* v2_30 - Dugas B.          - output real*8 cartesian coordinates
*
*object
*     See above ID
*
*arguments
*  Name        I/O                 Description
*----------------------------------------------------------------
* F_xyz_8      O    - coordinates in cartesian space
* F_lon        I    - longitudes of the grid in spherical coordinates
* F_lat        I    - latitudes of the grid in spherical coordinates
*
**
      integer i,j,k
      real*8, parameter :: ONE = 1.0d0, PI = 180.0d0
      real*8 deg2rad_8
*
*     __________________________________________________________________
*
      deg2rad_8 = acos(-ONE)/PI
*
      k=0
      do j=1,nj
      do i=1,ni
         k=k+1
         F_xyz_8(1,K) = cos(deg2rad_8*F_lat(j))*cos(deg2rad_8*F_lon(i))
         F_xyz_8(2,K) = cos(deg2rad_8*F_lat(j))*sin(deg2rad_8*F_lon(i))
         F_xyz_8(3,K) = sin(deg2rad_8*F_lat(j))
      enddo
      enddo
*
*     __________________________________________________________________
*
      return
      end
