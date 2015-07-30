      subroutine gc2rc(VN,vo,nn,AA)
!*********************************************************************
! author:	Xingxiu Deng
! date:		February 2007
! Object:	geographic cartesian coord to rotated cartesian coord 
!               matrix AA times a vector vo and get a new vector VN
!*********************************************************************
! input:        
! output:       
!*********************************************************************
      implicit none
!* input
      integer nn
      real*8 vo(nn), AA(nn,nn)

!* output
      real*8 VN(nn)

!* misc
      integer i,j
      logical :: debug=.false.

      if (debug) then
         do i=1,3
            print *, (AA(i,j),j=1,3)
         enddo
      endif

      DO i=1,nn
       VN(i)=0.0
        DO j=1,nn
         VN(i)=VN(i)+AA(i,j)*vo(j)
        ENDDO
      ENDDO

      if (debug) then
          print*, 'vo= ', (vo(i),i=1,nn)
          print*, 'VN= ', (VN(i),i=1,nn)
      endif

      RETURN
      END
