      SUBROUTINE sub_rotation(gdgem,lat,lon,ax,ay,AA,nx,ny,fstZfile) 
      implicit none
!*********************************************************************
! author: Sylvie Gravel, 30 March 2006
! changed to be a subroutine for its use in sharing factor generator
!                 --  Xingxiu Deng, April 2006
!*********************************************************************
! input:  fstZfile: FST file on Z grid
!         nx: grid dimension in x-direction
!         ny: grid dimension in y-direction
! output: gdgem: grid id
!         lat: latitude in geographical polar coord
!         lon: longitude in geographical polar coord
!         ax: longitude in rotated coord
!         ay: latitude in rotated coord
!*********************************************************************
      integer ezqkdef,gdgaxes,gdll
      integer gdgem, ier, nrecs
      integer fstinf, fstprm, fstouv, fnom, fstfrm, fclos
      external fstinf, fstprm, fstouv, fnom, fstfrm, fclos
      integer iun, ni,nj,nk,ip1,ip2,ip3,deet,npas,nbits,datyp,dateo
      integer ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,ex1,ex2,ex3
      integer datev, key
      character*2 nomvar
      character*1 typvar,grtyp
      character*8 etiket
      real pi
      real X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,XA(3),KA(3),KB(3)
      real alpha,c1,c2,c3,c4,beta
!
      integer ipos,i0,in,im,j0,jn,jm, i,j
! 
      integer nx, ny
      character (len=100) fstZfile

      real ax(nx), ay(ny)
      real lat(nx,ny), lon(nx,ny)
      real*8  AA(3,3), BB(2,3)

      pi = acos(-1.)

      iun = 1
      ier = fnom(iun,trim(fstZfile),'STD+RND',0)
      nrecs = fstouv(iun,'RND')

      datev  = -1
      etiket = '        '
      typvar = ' '
      nomvar = 'ME'
      ip1 = -1
      ip2 = -1
      ip3 = -1
!**   read fld info and data
      key=FSTINF(iun,NI,NJ,NK,datev,etiket,ip1,ip2,ip3,typvar,nomvar) 
      print*, 'NI=, NJ= ', NI, NJ
      print*, 'nx=, ny= ', nx, ny

      IF (NI .NE. nx .OR. NJ .NE.ny ) THEN
       print*, 'incorrect grid definition --- abort'
       print*, 'nx=, ny= ', nx, ny
       print*, 'NI=, NJ= ', NI, NJ
      ier = fstfrm(iun)
      ier = fclos(iun)
       stop
      ENDIF

      ier=FSTPRM(key,DATEO,DEET,NPAS,NI,NJ,NK,NBITS,DATYP,IP1, &
          IP2,IP3,TYPVAR,NOMVAR,ETIKET,GRTYP,IG1,IG2,IG3,      &
          IG4,SWA,LNG,DLTF,UBC,EX1,EX2,EX3)

!**   Define input grid         
      gdgem = ezqkdef(ni,nj ,grtyp,ig1, ig2, ig3, ig4, iun)

!**   Gets rotated latlon values for gdid "gdgem"
      ier = gdgaxes(gdgem, ax, ay)

!**   Gets geographical latlon values for gdid "gdgem"
      ier = gdll(gdgem, lat, lon)

      i0 = 1
      in = ni
      j0 = 1
      jn = nj
      print *, ' ax(i0),ay(j0) = ',ax(i0),ay(j0)
      print *, ' ax(in),ay(jn) = ',ax(in),ay(jn)
      print *, ' lat(i0,j0),lon(i0,j0) = ',lat(i0,j0),lon(i0,j0)
      print *, ' lat(in,jn),lon(in,jn) = ',lat(in,jn),lon(in,jn)

!**     cartesian coordinate of point (i0,j0) in rotated frame
!**
      XA(1) = cos(ay(i0)*pi/180.)*cos(ax(j0)*pi/180.)
      XA(2) = cos(ay(i0)*pi/180.)*sin(ax(j0)*pi/180.)
      XA(3) = sin(ay(i0)*pi/180.)
!*
!*     cartesian coordinate on geographical grid of following:
!*     AA(1) & AA(2) are on the longitude i0 of the rotated grid
!*     BB(1) & BB(2) are on the longitude j0 of the rotated grid
!*
      AA(1,1) = cos(lat(i0,j0)*pi/180.)*cos(lon(i0,j0)*pi/180.)
      AA(1,2) = cos(lat(i0,j0)*pi/180.)*sin(lon(i0,j0)*pi/180.)
      AA(1,3) = sin(lat(i0,j0)*pi/180.)
      AA(2,1) = cos(lat(i0,jn)*pi/180.)*cos(lon(i0,jn)*pi/180.)
      AA(2,2) = cos(lat(i0,jn)*pi/180.)*sin(lon(i0,jn)*pi/180.)
      AA(2,3) = sin(lat(i0,jn)*pi/180.)
      BB(1,1) = cos(lat(in,j0)*pi/180.)*cos(lon(in,j0)*pi/180.)
      BB(1,2) = cos(lat(in,j0)*pi/180.)*sin(lon(in,j0)*pi/180.)
      BB(1,3) = sin(lat(in,j0)*pi/180.)
      BB(2,1) = cos(lat(in,jn)*pi/180.)*cos(lon(in,jn)*pi/180.)
      BB(2,2) = cos(lat(in,jn)*pi/180.)*sin(lon(in,jn)*pi/180.)
      BB(2,3) = sin(lat(in,jn)*pi/180.)
!*
!*     the 2 intersections of the above longitudes are the numerical
!*     poles of the rotated grid.
!*
!*     KA is the cartesian coordinate of the normal to the great circle
!*     that passes through both AA(1) and AA(2)

      KA(1) = AA(1,2)*AA(2,3)-AA(1,3)*AA(2,2)
      KA(2) = AA(1,3)*AA(2,1)-AA(1,1)*AA(2,3)
      KA(3) = AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1)
      c1 = KA(1)*KA(1)+KA(2)*KA(2)+KA(3)*KA(3)
      KA(1) = KA(1)/c1
      KA(2) = KA(2)/c1
      KA(3) = KA(3)/c1
!*     KB is the cartesian coordinate of the normal to the great circle
!*     that passes through both BB(1) and BB(2)

      KB(1) = BB(1,2)*BB(2,3)-BB(1,3)*BB(2,2)
      KB(2) = BB(1,3)*BB(2,1)-BB(1,1)*BB(2,3)
      KB(3) = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
      c1 = KB(1)*KB(1)+KB(2)*KB(2)+KB(3)*KB(3)
      KB(1) = KB(1)/c1
      KB(2) = KB(2)/c1
      KB(3) = KB(3)/c1
!*
!*     by definition:
!*     the numerical pole X0,Y0,Z0 is a point that is perpendicular to both 
!*     KA, and KB since it belong to both great circle, the one defined
!*     by AA(1)& AA(2), and the one defined by BB(1) & BB(2)
!*     this relationship allows us to obtain X0,Y0,Z0 the numerical pole
!*
      c1 = -KA(2)/KA(1)
      c2 = -KA(3)/KA(1)
      c3 = -(c2*KB(1)+KB(3))/(c1*KB(1)+KB(2))
      c4 = (c1*c3 +c2)
      Z0 = 1.+c3*c3+c4*c4
      Z0 = sqrt(1./Z0)
      X0 = c4*Z0
      Y0 = c3*Z0
!*
!*     we define 2 vectors r1 and r2 perpendicular to each other
!*     and in the equatorial plane having X0,Y0,Z0 for a normal
!*
      c2 = -Z0/X0
      Z1 = 1.+c2*c2
      Z1 = sqrt(1./Z1)
      Y1 = 0.   ! arbitrary choice
      X1 = c2*Z1
!*
      c1 = -Y0/X0
      c2 = -Z0/X0
      c3 = -(X1*c2+Z1)/(X1*c1)
      c4 = c1*c3+c2
      Z2 = c4*c4+c3*c3+1.
      Z2 = sqrt(1./Z2)
      Y2 = c3*Z2
      X2 = c4*Z2
!*
!*     check to see if cross product of r1 and r2 = (X0,Y0,Z0)
!*     if necessary, change the sign of r2
!*
      if ( Z1*Y2*X0 .gt. 0.) then
         X2 = -X2
         Y2 = -Y2
         Z2 = -Z2
      endif
!*    
!*     the center of the domain, whose cartesian coordinates are the
!*     top row of the rotation matrix must be a linear combinaison
!*     of r1 and r2 since r1 and r2 are in the equatorial plane of
!*     the rotated grid.  We will store this vector in AA(1,i)
!*
!*     similarly the second row of the matrix must also be a combinaison
!*     of r1 and r2, and chosen such that the dot product of AA(1,i) and 
!*     AA(2,i) is zero, and their cross product = X0,Y0,Z0
!*
!*     using these constraints and the fact that we know the
!*     coordinates of the point i0,j0 in both the rotated and
!*     geographical system, we can compute the values of the 1st and
!*     2nd row of the rotational matrix.  The last row is given by
!*     X0, Y0, Z0.
!*
      c1 = (AA(1,1)*X1 + AA(1,2)*Y1 + AA(1,3)*Z1)
      c2 = (AA(1,1)*X2 + AA(1,2)*Y2 + AA(1,3)*Z2)

      alpha = (c2*XA(1)-c1*XA(2))/(c1*c1+c2*c2)
      beta = (XA(1)-c2*alpha)/c1

      print *, ' rotation matrix '
      AA(1,1) = beta*X1 + alpha*X2
      AA(1,2) = beta*Y1 + alpha*Y2
      AA(1,3) = beta*Z1 + alpha*Z2
      AA(2,1) = alpha*X1 - beta*X2
      AA(2,2) = alpha*Y1 - beta*Y2
      AA(2,3) = alpha*Z1 - beta*Z2
      AA(3,1) = X0
      AA(3,2) = Y0
      AA(3,3) = Z0
      if ( (AA(1,3)*AA(2,2)-AA(1,2)*AA(2,3))*AA(3,1) .gt. 0.) then
         AA(2,1) = -AA(2,1)
         AA(2,2) = -AA(2,2)
         AA(2,3) = -AA(2,3)
      endif

      do i=1,3
      print*, (AA(i,j),j=1,3)
      enddo

! close fst file
      ier = fstfrm(iun)
      ier = fclos(iun)

      RETURN
      END SUBROUTINE sub_rotation
