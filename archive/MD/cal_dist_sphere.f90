      subroutine cal_dist_sphere(xyz1_8,xyz2_8,r)
!**********************************************************
!* author: Xingxiu Deng, Mar., 2006
!* purpose: to calculate the spherical distance between two
!* points at the earth surface (defined by Cartesian coordinates)
!**********************************************************
      implicit none
!*input 
      real*8 xyz1_8(3),xyz2_8(3)

!*output (r in m), can't be in real*8, since the calling subroutine
!     use real only
      real r

!*misc, angrad (in radiance), rearth (the radius of the Earth in km)
      real*8 rs,angrad,rearth

      rearth=6378.135

      rs=sqrt((xyz2_8(1)-xyz1_8(1))*(xyz2_8(1)-xyz1_8(1))+ &
              (xyz2_8(2)-xyz1_8(2))*(xyz2_8(2)-xyz1_8(2))+ &
              (xyz2_8(3)-xyz1_8(3))*(xyz2_8(3)-xyz1_8(3)))
      angrad=2.*asin(rs*0.5)
      r=rearth*angrad*1000.

      RETURN
      END
