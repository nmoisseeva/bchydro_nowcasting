      SUBROUTINE MD_sph(n_HA,name_HA,HAlonrot,HAlatrot,HAxyz_rot_8,z_HA, &
                             ni,nj,lonrot,latrot,xyzgr_rot_8,Zdaughter,  &
                             Rxy,dx,dy,a,b,zref1,zref2,avehpbl,          &
                             Sdaughter_3D,rg_3D,i_HAx,j_HAy,             &
                             Rnearest_HA,Znearest_HA,Sdaughter_3Dp)
!*********************************************************************
! author:	Xingxiu Deng
! date:		2003 - 2004 (original version)
!               January - April 2006 (major modification for spheric coords)
!               March 2007 (documentation, major clean up and adding
!                  a new subroutine to shorten the main subroutine)
! Object:	to calculate sharing factors for the entire domain or
!               within an influencial region from each source point using 
!               the mother-daughter (MD) approach
!*********************************************************************
! Algorithm:
!               For detail information about the MD approach, please
!               refer to the following papers:
!
!    Deng, X. and R. Stull, 2005: A mesoscale analysis method for surface
!    potential temperature in mountainous and coastal terrain. Mon. Wea.
!    Rev., 133, 389-408.
!
!    Deng, X. and R. Stull, 2007: Assimilating Surface Weather Observations
!    from Complex Terrain into a High-Resolution Numerical Weather Prediction
!    Model. Mon. Wea. Rev., 135, 1037-1054.
!*********************************************************************
! input:        
!       for source points
!               n_HA: No. of source points
!               name_HA: name or character id of the source points
!               HAlatrot: rotated latitude of the source points
!               HAlonrot: rotated longitude of the source points
!               HAxyz_rot_8: rotated Cartesian coords(3) of the source points
!               z_HA: elevation of the source points (m)
!
!       for working grid: 
!               ni,nj: dimension of the working grid
!               lonrot: rotated longtitude of the working grid
!               latrot: rotated latitude of the working grid
!               xyzgr_rot_8: rotated Cartesian coords(3) of the working grid
!               Zdaughter: topography of the working grid (m)
!
!               Rxy: Radius of influence (km)
!               dx: grid spacing in x-dir
!               dy: grid spacing in y-dir
!       derived: R_gaussian: same as Rxy but in meters
!
!       for MD free parameters:
!               a, b: parameters controling analysis decorrelation rate
!               zref1: terrain-following BL-depth parameter (m)
!               zref2: level-top BL-depth parameter (m) 
!               avehpbl: approximated average BL depth (m)
!
! output:       Sdaughter_3D: sharing factors at grid points from all 
!                             source points/stations
!               Sdaughter_3Dp: same as Sdaughter_3D but after Gaussian dropoff 
!               rg_3D: circuitous travel distance (CTD) from source point to 
!                      the GP (m) (=ini_rg=99999.E03 if the source point has no
!                       influence on the GP)
!       Following output are used for obtaining sharing factor from one
!       source point to another source point rather than to a GP
!               i_HAx: i-index of the source point's nearest GP
!               j_HAy: j-index of the source point's nearest GP
!               Rnearest_HA: distance between the GP at (i_HAx, j_HAy) and 
!                            the source point
!               Znearest_HA: height diff between the GP at (i_HAx, j_HAy) and
!                            the source point
!*********************************************************************
! Modification history:
! --- added ibegin,iend,jbegin and jend to define the true range of 
!     grid points for updating Sdaughter, in order to reduce CPU time 
!     for all iterations (May. 2003)
! --- added refinement for coastal terrain (June. 2003)
! --- added influence region (R_inf) to limit ibegin,iend,jbegin and jend,
!     again in order to reduce CPU time. The influence region should
!     be larger than the influence region used in data analysis
!     process (June, 2004)
! --- added refinement for mountain-top observations (July, 2004)
! --- removed lat/lon (HAlat/HAlon, latgr/longr) for source points
!     and other grid points, remove PI_nps.... ( Jan. 23, 2006)
! --- added dy in order to work for different dx and dy (Jan. 24, 2006)
! --- added if_limit to check whether or not apply influence region when
!     calculating sharing factors (false: no; true: yes) (Jan. 24, 2006)
! --- removed the two refinements temporarily     (Jan. 26, 2006)
! --- applied to the GEM rotated lat-lon coordinate, HAx --> HAlonrot
!     HAy --> HAlatrot, xgr, ygr --> lonrot, latrot. added HAxyz_rot_8
!     and xyzgr_rot_8 for the use in calculating distance. (Apr., 2006)
!     xyzgr_8 and HAxyz_8 can now only be used to calculate distance, and 
!     HAlatrot/HAlonrot, and latrot/lonrot have to be used for comparing 
!     relative position of grid points (Apr. 2006).
! --- Document input and output variables, major clean up, improved
!     overall document on the code, and shortened the main subroutine
!     by adding a new subroutine check_mother (Mar.- Apr. 2007)
! Note: if_limit = .false. for 3D-var, if = .true., need to add 
!      calculation of R_inf (Mar. 2007)
!**********************************************************
         IMPLICIT none

!*** input variables:
!--------------------------------------------------
!* source points/stations
         INTEGER :: n_HA
         CHARACTER (len=5) name_HA(n_HA)
         REAL :: HAlonrot(n_HA), HAlatrot(n_HA)
         REAL*8 :: HAxyz_rot_8(n_HA,3)
         REAL :: z_HA(n_HA)
!* derived input for source points/stations
         REAL :: HAx(n_HA),HAy(n_HA)
         REAL :: S_HA(n_HA)         ! always = 1.

!* working grid
         INTEGER :: ni,nj
         REAL :: lonrot(ni),latrot(nj) 
         REAL*8 :: xyzgr_rot_8(ni,nj,3)
         REAL Zdaughter(ni,nj)
         REAL Rxy
         REAL dx, dy

!* derived input for working grid
         REAL :: xgr(ni,nj),ygr(ni,nj)
         REAL :: R_gaussian

!* free parameters
         REAL :: a,b
         REAL :: zref1, zref2
         REAL :: avehpbl
             
!*** output variables:
!--------------------------------------------------
         REAL Sdaughter_3D(ni,nj,n_HA),Sdaughter_3Dp(ni,nj,n_HA)
         REAL rg_3D(ni,nj,n_HA)
         INTEGER :: i_HAx(n_HA)
         INTEGER :: j_HAy(n_HA)
         REAL :: Rnearest_HA(n_HA)
         REAL :: Znearest_HA(n_HA)

         INTEGER ncycl_3D(ni,nj,n_HA)

!*** intermediate variables:
!--------------------------------------------------
         INTEGER :: index_hax     ! x-dir index of the source point's nearest GP
         INTEGER :: index_hay     ! y-dir index of the source point's nearest GP
         REAL :: rmin             ! distance between the source point and its nearest GP
         REAL :: S_mother         ! sharing factor for a temporary mother
         REAL :: Z_mother         ! terrain height for a temporary mother
         INTEGER :: ipt, jpt      ! index of the GP to the south west of the source point
         INTEGER :: ireturn       ! indicator (=0, the source point is inside the grid domain;
                                  !            <0, the source point is outside the grid domain)
         INTEGER :: icycle        ! counter of iterative cycle
         INTEGER :: ibegin        ! true range of GPs to update Sdaughter.
         INTEGER :: iend          ! they depend on icycle.
         INTEGER :: jbegin        ! ....
         INTEGER :: jend          ! Added to reduce CPU time for all iterations.
         INTEGER :: ibegin_min, iend_max
         INTEGER :: jbegin_min, jend_max
         REAL :: tmp1,tmp2        ! temporary array
         REAL :: tmp3,tmp4        ! temporary array
         REAL R_inf               ! Radius of influence region(m), used to limit ibegin,iend,jbegin,jend.
         PARAMETER (R_inf=300.E03) ! if Rxy=90.km, cut-off radii in ADAS is 273.14km. then R_inf=300 km
         LOGICAL :: if_limit = .false. ! logical variable, to tell whether or not apply influence region
                                       ! when calculating sharing factor (true: yes; false:no)

         INTEGER :: indm           ! order (or sequence) of the mother around the grid point
         REAL :: s(9)             ! sharing factor from every possible mothers for a GP
         REAL :: r(9)             ! distance from every possible mothers to a GP
         REAL :: imother(9)       ! i index for every possible mothers
         REAL :: jmother(9)       ! j index for every possible mothers 

         INTEGER ncycl(ni,nj),ncycl0(ni,nj)  ! variables to memorize no. of cycle for each GP 
         REAL rg(ni,nj)                      ! minimum travel distance from source point to the GP (m)
         REAL rg0(ni,nj)                     ! minimum travel distance from SP to a possible mother
         REAL rg_out(ni,nj),rg_out0(ni,nj)   ! similar to rg and rg0, but initialized to be 
                                             ! 99999.E03(ini_rg) rather than 0., rg_3D = rg_out
         REAL r_min(ni,nj)                   ! minimum travel distance from the source point to GP
                                             ! minimum of r_smax if there are several mother that has 
                                             ! the same largest contributions 
         REAL r_smax(ni,nj,9)                ! the distance from the source point to the GP via the
                                             ! mother point that has the largest contributions (m)
         REAL :: ini_rg                      ! m
         PARAMETER (ini_rg = 99999.E03)     

         REAL Smother(ni,nj)      ! sharing factor for every possible mother
         REAL Sdaughter(ni,nj)    ! sharing factor for daughters
         REAL Stmp(ni,nj)         ! temporary array for Sdaughter for comparison in iterative process.

!*** Misc variables
!--------------------------------------------------
         INTEGER :: ista, i,j, ii,jj,i1,j1,i2,j2,i0,j0
         INTEGER :: k,kk
         LOGICAL :: debug = .true.
         REAL :: xgrmin,xgrmax,ygrmin,ygrmax
!--------------------------------------------------

         print*, '---- EXECUTE MD_sph -------------'

!*** switch new variables to old ones
          DO ista = 1, n_HA
           HAx(ista) = HAlonrot(ista)
           HAy(ista) = HAlatrot(ista)
          ENDDO
          DO j=1,nj
          DO i=1,ni
           xgr(i,j)=lonrot(i)
           ygr(i,j)=latrot(j)
          ENDDO
          ENDDO
          R_gaussian = Rxy * 1.0E03

!*** Initialize arrays
          Sdaughter_3D = 0.
          Sdaughter_3Dp = 0.
          rg_3D = ini_rg
          ncycl_3D = 0
          Rnearest_HA = 0.
          Znearest_HA = 0.
          i_HAx = 0
          j_HAy = 0
          S_HA = 1.

!*** find xgrmin,xgrmax,ygrmin,ygrmax
         xgrmin = 1.0E+8
         xgrmax = -1.0E+8
         ygrmin = 1.0E+8
         ygrmax = -1.0E+8

         DO j =1,nj
           IF (xgr(1,j) .LT. xgrmin) xgrmin =  xgr(1,j)
           IF (xgr(ni,j) .GT. xgrmax) xgrmax =  xgr(ni,j)
         ENDDO      
         DO i =1,ni
           IF (ygr(i,1) .LT. ygrmin) ygrmin =  ygr(i,1)
           IF (ygr(i,nj) .GT. ygrmax) ygrmax =  ygr(i,nj)
         ENDDO      

!*** loop over all the source points
      DO 1000 ista = 1, n_HA

!*** Initialize the following two arrays for each source point
        Smother = 0.
        Sdaughter = 0.

!*** find the GP to the south and west of the source point
        IF ( HAx(ista) .LT. xgrmin .OR. HAx(ista) .GT. xgrmax .OR.    &
             HAx(ista) .LT. xgrmin .OR. HAx(ista) .GT. xgrmax  ) THEN
          ireturn = -1
          if (debug) then
           print*,'xgrmin,xgrmax,ygrmin,ygrmax',xgrmin,xgrmax,ygrmin,ygrmax
           print*,'HAx,HAy',HAx(ista),HAy(ista)
          endif

          print*, name_HA(ista), ' is outside the domain'
          Sdaughter_3D = -999.
          Sdaughter_3Dp = -999.
          rg_3D = -999.
          ncycl_3D = -999
          Rnearest_HA = -999.
          Znearest_HA = -999.
          i_HAx = -999
          j_HAy = -999
          GO TO 1000
        ELSE
           CALL findlc(ni,nj,xgr,ygr,HAx(ista),HAy(ista),ipt,jpt,ireturn)
           if (debug) then
            print*, 'ipt,jpt,ireturn',ipt,jpt,ireturn
           endif
        ENDIF

!******************************************************************
! NOTE (Xingxiu Deng, 2003):
!  When the source point is not collocated with a grid point, 
!  an approximation is applied (the source point is approximated
!  to be collocated with the grid point that is nearest to the
!  source point in the sense of minimum elevation difference. 
!  Need to account for elevation difference between the source point
!  and its nearest GP. So the sharing factor to the nearest GP is reduced  
!  by elevation difference between them using a Gaussian drop-off with
!  a standard deviation=avehpbl. And then Z_HA should be corrected 
!  to Zdaughter(index_hax,index_hay))
!******************************************************************

!*** Find nearest GP to the source point if it is not collcated with a GP
        IF (HAx(ista) .NE. xgr(ipt,jpt) .AND.  &
            HAy(ista) .NE. ygr(ipt,jpt)) THEN
          CALL find_index_HA_sph(ipt,jpt,HAxyz_rot_8(ista,:),Z_HA(ista), &
                              ni,nj,xyzgr_rot_8,Zdaughter, &
                              index_HAx,index_HAy,rmin)
          i_HAx(ista) = index_HAx
          j_HAy(ista) = index_HAy
          Rnearest_HA(ista) = rmin
          Znearest_HA(ista) = Z_HA(ista)-Zdaughter(index_hax,index_hay)
          Smother(index_hax,index_hay) = S_HA(ista)*         &
                  exp(-0.5*Znearest_HA(ista)**2/avehpbl**2)
          Sdaughter(index_hax,index_hay) = Smother(index_hax,index_hay)
          Z_HA(ista) = Zdaughter(index_hax,index_hay)
        ELSE
          index_HAx=ipt
          index_HAy=jpt
          i_HAx(ista) = index_HAx
          j_HAy(ista) = index_HAy
          rmin=0.
          Rnearest_HA(ista) = rmin
          Znearest_HA(ista) = 0.
          Smother(index_HAx,index_HAy) = S_HA(ista)
          Sdaughter(index_HAx,index_HAy) = Smother(index_hax,index_hay)
!          Z_HA(ista) = Z_HA(ista)
        ENDIF

!*** Memorize sharing factors at the begining for comparison with those 
!    after one iteration
        Stmp = Sdaughter

!*** 1st cycle 
        icycle = 1

!*** Initialize ncycl and ncycl0 for every grid point
        ncycl = 0
        ncycl0 = ncycl
!*** Initialize rg and rg0 for every grid point
        rg = 0.
        rg0 = rg
        rg_out = ini_rg
        rg_out(index_hax,index_hay) = 0.
        rg_out0 = rg_out

!*** initialize ibegin, iend, jbegin, jend
        ibegin_min = index_HAx 
        iend_max = index_HAx 
        jbegin_min = index_HAy 
        jend_max = index_HAy 

        ibegin = ibegin_min
        iend = iend_max
        jbegin = jbegin_min
        jend = jend_max

10      CONTINUE
!*************************************************************************
! initialize s and r for every cycle
!*************************************************************************
        s = 0.
        r = 0.

!*************************************************************************
! loop over ibegin, iend, jbegin, jend
!*************************************************************************
        IF (if_limit) THEN
         tmp1=sqrt(((ibegin_min-index_HAx)*dx)**2)   
         tmp2=sqrt(((jbegin_min-index_HAy)*dy)**2)   
         tmp3=sqrt(((iend_max-index_HAx)*dx)**2)    
         tmp4=sqrt(((jend_max-index_HAy)*dy)**2)   
         IF (tmp1 .LT. R_inf   .OR.  tmp2 .LT. R_inf   .OR.   &
             tmp3 .LT. R_inf   .OR.  tmp4 .LT. R_inf ) THEN    
           ibegin_min = index_HAx - icycle
           iend_max = index_HAx + icycle
           jbegin_min = index_HAy - icycle
           jend_max = index_HAy + icycle
         ENDIF
        ELSE
         ibegin_min = index_HAx - icycle
         iend_max = index_HAx + icycle
         jbegin_min = index_HAy - icycle
         jend_max = index_HAy + icycle
        ENDIF
        ibegin = ibegin_min
        iend = iend_max
        jbegin = jbegin_min
        jend = jend_max

        IF ( ibegin .LT. 1 ) ibegin = 1
        IF ( iend .GT. ni ) iend = ni
        IF ( jbegin .LT. 1 ) jbegin = 1
        IF ( jend .GT. nj ) jend = nj

        DO j = jbegin, jend
        DO i = ibegin, iend

! loop over all the mothers around each daughter         
         DO indm = 1, 8
          IF ( indm .eq. 1) THEN  
           i1 = i+1
           j1 = j+1
          ELSEIF ( indm .eq. 2) THEN  
           i1 = i+1
           j1 = j
          ELSEIF ( indm .eq. 3) THEN  
           i1 = i+1
           j1 = j-1
          ELSEIF ( indm .eq. 4) THEN  
           i1 = i
           j1 = j-1
          ELSEIF ( indm .eq. 5) THEN  
           i1 = i-1
           j1 = j-1
          ELSEIF ( indm .eq. 6) THEN  
           i1 = i-1
           j1 = j
          ELSEIF ( indm .eq. 7) THEN  
           i1 = i-1
           j1 = j+1
          ELSEIF ( indm .eq. 8) THEN  
           i1 = i
           j1 = j+1
          ENDIF
          call  check_mother(ista,z_HA,n_HA,ni,nj,Smother,Zdaughter, &
                            xyzgr_rot_8,i1,j1,i,j,a,b,zref1,zref2,indm, &
                            imother,jmother,s,r)
         ENDDO

         indm = 9              
         s(indm) = Sdaughter(i,j)
         r(indm) = 0.
         imother(indm) = i
         jmother(indm) = j

!*** Keep the strongest daughter at the given point (i,j)        
         Sdaughter(i,j) = amax1(s(9),s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8))

!*** Find the nearest mother that has maximum contribution to the grid point
         IF (Sdaughter(i,j) .NE. 0. ) THEN
          kk = 0
          DO k =1, 9
           IF (Sdaughter(i,j) .EQ. s(k)) THEN
            kk = kk + 1
            i2 = imother(k)
            j2 = jmother(k)
            r_smax (i,j,kk) = r(k)+ rg0(i2,j2)
           ENDIF
          ENDDO
          r_min(i,j)= r_smax(i,j,1)
          DO k = 1,kk
           IF( r_smax(i,j,k) .LT. r_min(i,j)) r_min(i,j) = r_smax(i,j,k)
          ENDDO
!        ELSE
!         no changes to travelling dictance
         ENDIF

        ENDDO    ! end of i
        ENDDO    ! end of j

!*************************************************************************
! At this point, one cycle is finished for every grid point (i,j)
! Following lines are for debuging and getting ready for next cycle
!*************************************************************************
        DO j = 1, nj
        DO i = 1, ni
          IF (abs(Sdaughter(i,j)-Stmp(i,j)) .GE. 0.01) THEN
           ncycl(i,j) = icycle 
           rg(i,j) =  r_min(i,j)
           rg_out(i,j) =  r_min(i,j)
          ELSE
           ncycl(i,j) = ncycl0(i,j)
           rg(i,j) = rg0(i,j)
           rg_out(i,j) = rg_out0(i,j)
          ENDIF
        ENDDO
        ENDDO

!*** Update Smother after one cycle
        Smother = Sdaughter

!*** Update ncycl0, and rg0
        ncycl0 = ncycl
        rg0 = rg
        rg_out0 = rg_out
        
!*** Check whether or not Sdaughter over entire domain update
!*** compared to the value (Stmp) of the previous cycle
        DO j =1, nj
        DO i =1, ni
          IF (abs(Sdaughter(i,j)-Stmp(i,j)) .GE. 0.01) GO TO 888
        ENDDO
        ENDDO

        GO TO 999

888     CONTINUE
        Stmp = Sdaughter
        icycle = icycle + 1 
        GO TO 10

!*** exit cycle
999     CONTINUE

        DO j =1,nj
        DO i =1,ni
          IF (abs(Zdaughter(i,j)-Z_HA(ista)).GE.Zref2) rg_out(i,j)= ini_rg ! m
        ENDDO
        ENDDO

        DO j =1,nj
        DO i =1,ni
          Sdaughter_3D(i,j,ista) = Sdaughter(i,j)
          ncycl_3D(i,j,ista) = ncycl(i,j)
          rg_3D(i,j,ista) = rg_out(i,j) + Rnearest_HA(ista)      ! meters
        ENDDO
        ENDDO

!*************************************************************************
! sharing factor after Gaussian drop-off (for plot):
! rg_out (meters), Rnearest_HA(meters), R_gaussian(meters)
!-------------------------------------------------------------------------
        DO j0 = 1, nj
        DO i0 = 1, ni
           Sdaughter_3Dp(i0,j0,ista) = Sdaughter(i0,j0) *     &
              exp(-0.5*((rg_out(i0,j0)+Rnearest_HA(ista))/R_gaussian)**2)
        ENDDO
        ENDDO
!*************************************************************************

 1000   CONTINUE       ! end of ista

        RETURN
        END

        SUBROUTINE check_mother(ista,z_HA,n_HA,ni,nj,Smother,Zdaughter, &
                               xyzgr_rot_8,i1,j1,i,j,a,b,zref1,zref2,indm, &
                                imother,jmother,s,r)
!***********************************************************
! author: Xingxiu Deng
! date: March 2007
! purpose: to check each mother's grid index and call SUBROUTINE 
!          mother and distance_sph to shorten the main SUBROUTINE
!          MD_sph
!***********************************************************
! input variables:
! ista: index of the source point
! z_HA: elevation of the source points (m)
! n_HA: No. of source points
! ni : x-dir dimension size
! nj : y-dir dimension size
! Smother : sharing factors of mothers
! Zdaughter: elevation of the daughter (m)
! xyzgr_rot_8: cartesian coord of grid points
! i1,j1: grid index of a mother
! i,j: grid index of a daughter
! a, b, zref1,zref2: MD's free parameters
! indm: order (or sequence) of the mother around the grid point
!
! output variable:
! imother:  i index for every possible mothers
! jmother:  j index for every possible mothers
! s :  sharing factors at the point(id,jd) from the mother (im,jm)
! r :  distance from the point(id,jd) to the mother (im,jm)
!-----------------------------------------------------------
         IMPLICIT NONE
         INTEGER :: n_HA, ista
         REAL :: z_HA(n_HA)

         INTEGER :: ni,nj,i1,j1,i,j,indm
         REAL :: Zdaughter(ni,nj),Smother(ni,nj)
         REAL*8 :: xyzgr_rot_8(ni,nj,3) 
         REAL :: a,b,zref1,zref2

         REAL :: imother(9), jmother(9)
         REAL :: s(9), r(9)
         REAL :: ss, rr

         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(indm) = i1
          jmother(indm) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,ss,  &
                      Zdaughter,i,j,a,b,zref1,zref2)
          s(indm) = ss
          CALL distance_sph(rr,i1,j1,i,j,xyzgr_rot_8,ni,nj)
          r(indm) = rr
         ELSE
          s(indm) = 0.
          r(indm) = 0.
          imother(indm) = i
          jmother(indm) = j
         ENDIF

         RETURN
         END SUBROUTINE check_mother

         SUBROUTINE mother(ista,z_HA,n_HA,ni,nj,Smother,im,jm,ss,   &
                             Zdaughter,id,jd,a,b,zref1,zref2)
!***********************************************************
! author: Xingxiu Deng
! date: Jan. 20, 2003
! purpose: to clculate sharing factors of a grid point (id,jd) 
!          from one of its mother (Smother(im,jm))
!***********************************************************
! input variables:
! ista: index of the source point
! z_HA: elevation of the source points (m)
! n_HA: No. of source points
! ni : x-dir dimension size
! nj : y-dir dimension size
! Smother : sharing factors of mothers
! im,jm: grid index of a mother
! Zdaughter: elevation of the daughter (m)
! id,jd: grid index of a daughter
! a, b, zref1,zref2: MD's free parameters
!
! output variable:
! ss :  sharing factors at the point(id,jd) from the mother (im,jm)
!
! intermediate:
! Z_mother : elevation of the mother (m)
!-----------------------------------------------------------
         IMPLICIT NONE
         INTEGER :: n_HA, ista
         REAL :: z_HA(n_HA)

         INTEGER :: ni,nj,im,jm,id,jd,ind
         REAL :: Zdaughter(ni,nj),Smother(ni,nj)
         REAL :: a,b,zref1,zref2

         REAL :: ss

         REAL :: Z_mother
      
! misc variables
         INTEGER :: i,j

         Z_mother = Zdaughter(im,jm)

         ss = Smother(im,jm) *                        &
          (1.-((abs(Z_mother - Zdaughter(id,jd)))/zref1)**a) *  &
          (1.-((abs(z_HA(ista) - Zdaughter(id,jd)))/zref2)**b)

         IF ( abs(z_HA(ista) - Zdaughter(id,jd)) .GT. zref2 ) then
           ss = 0.0
         ENDIF
         IF ( abs(Z_mother - Zdaughter(id,jd)) .GT. zref1 ) then
           ss = 0.0
         ENDIF
 
         RETURN
         END

         SUBROUTINE distance_sph(rr,im,jm,id,jd,xyzg_8,ni,nj)
!**********************************************************
! author: Xingxiu Deng
! date: April, 2006
! purpose: to calculate distance of a grid point (id,jd) 
!          to one of its mother (im,jm) at the surface of the sphere
! Based on SUBROUTINE distance (Xingxiu Deng, Jan. 27, 2003) but
! modified for GEM grid
!**********************************************************
         IMPLICIT NONE
!*** input
         INTEGER :: im,jm       ! grid index of a mother
         INTEGER :: id,jd       ! grid index of a daughter
         INTEGER :: ind          ! order (or sequence) of the mother around the GP
!*** output
         REAL :: rr           ! distance from the point(id,jd) to the mother (im,jm)
      
!*** intermediate
         REAL*8 :: xyzd_8(3)    ! cartesian coord for the grid point (m)
         REAL*8 :: xyzm_8(3)    ! cartesian coord for the mother (m)
         INTEGER :: ni,nj       ! grid dimension size
         REAL*8 :: xyzg_8(ni,nj,3) ! cartesian coord of grid points
         
         integer :: k
         
         DO k = 1, 3
          xyzd_8(k) = xyzg_8(id,jd,k)
          xyzm_8(k) = xyzg_8(im,jm,k)
         ENDDO
         
         call cal_dist_sphere(xyzd_8, xyzm_8,rr)

         RETURN
         END

        SUBROUTINE find_index_HA_sph(ipt,jpt,xyzpt_8,zpt,nx,ny,   &
                                  xyzg_8,zg,index_xpt,index_ypt,rmin)
!**********************************************************
! author: Xingxiu Deng, April 2006
! purpose: to find a nearest GP to the source point (xpt,ypt)
!          given the GP (ipt, jpt) to the south and west of (xpt, ypt)
!
!       - based on find_index_HA (Xingxiu Deng, Apr., 2003) which finds
!         the nearest GP in terms of minimum elevation difference. rmin
!         is defined as the horizontal distance between the source point
!         and its nearest GP. 
!
!       - NOW rmin is defined for GEM grid as great circle distance 
!         between the source point and its nearest GP.
!**********************************************************
        IMPLICIT NONE
!*** input         
        INTEGER :: ipt     ! i index to the west of desired
                           ! location (xpt, ypt)
        INTEGER :: jpt     ! j index to the south of desired
                           ! location (xpt,ypt)
        REAL :: zpt        ! height of desired location
        REAL*8 :: xyzpt_8(3) ! rotated cartesian coords of desired location
         
        INTEGER :: nx,ny   ! Grid dimensions 
        REAL :: zg(nx,ny)  ! z coordinate of grid points in 
                           ! physical space (m)
        REAL*8 :: xyzg_8(nx,ny,3)  ! rotated cartesian coords of grid points

        INTEGER :: index_xpt  ! i index that is closest to the desired
                              ! location (xpt, ypt)
        INTEGER :: index_ypt  ! j index that is closest to the desired
                              ! location (xpt, ypt)
!*** output
        REAL :: rmin          ! great circle distance between (xpt, ypt) and its nearest GP.

!*** Misc
        REAL :: zr(4)           
        REAL :: zrmin
        INTEGER :: k, kz

!*---------- Begin executable code ----------------------
         zr(1) = SQRT ( (zg(ipt,jpt) - zpt)**2 )
         zr(2) = SQRT ( (zg(ipt,jpt+1) - zpt)**2)
         zr(3) = SQRT ( (zg(ipt+1,jpt+1) - zpt)**2 )
         zr(4) = SQRT ( (zg(ipt+1,jpt) - zpt)**2)

         zrmin = amin1 ( zr(1), zr(2), zr(3),zr(4) )
         DO k = 1, 4
           IF ( zrmin .eq. zr(k) ) kz = k
         ENDDO
         IF ( kz .eq. 1 ) THEN
          index_xpt = ipt
          index_ypt = jpt
         ELSE IF ( kz .eq. 2 ) THEN
          index_xpt = ipt
          index_ypt = jpt + 1
         ELSE IF ( kz .eq. 3 ) THEN
          index_xpt = ipt + 1
          index_ypt = jpt + 1
         ELSE IF ( kz .eq. 4 ) THEN
          index_xpt = ipt + 1
          index_ypt = jpt 
         ENDIF
         call cal_dist_sphere(xyzg_8(index_xpt,index_ypt,:),xyzpt_8,rmin)

         print*, 'rmin,index_xpt,index_ypt',rmin,index_xpt,index_ypt
         print*, 'zg', zg(index_xpt,index_ypt)

         RETURN
         END SUBROUTINE find_index_HA_sph
