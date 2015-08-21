      subroutine rc2rp(latrot,lonrot,xyz_rot_8,nn)
!*********************************************************************
! author:	Xingxiu Deng
! date:		February 2007
! Object:	rotated cartesian coord to rotated polar coord
!*********************************************************************
! input:        
! output:       
!*********************************************************************
      implicit none
!* input
      integer nn
      real*8 xyz_rot_8(nn)

!* output
      real latrot, lonrot

!* misc
      real*8 pi
      real*8 rad2deg
      real my_aa

      pi = acos(-1.)
      rad2deg = 180./pi

      lonrot = rad2deg * (atan(xyz_rot_8(2)/xyz_rot_8(1)))
      latrot= rad2deg * (pi/2. - acos(xyz_rot_8(3)))

!* here lonrot is always smaller than 180., but lonrot from sub_rotation
!* might be bigger than 180. How to solve this problem?
      if (xyz_rot_8(1) .LT. 0.) then
         my_aa = 180.
      else if (xyz_rot_8(2) .GT. 0.) then
         my_aa = 0.
      else if (xyz_rot_8(2) .LT. 0.) then
         my_aa = 360.
      endif
      lonrot = my_aa + lonrot

      RETURN
      END
