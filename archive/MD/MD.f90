      SUBROUTINE cal_S_real_3D(n_HA,name_HA,HAx,HAy,z_HA,      &
                             ni,nj,xgr,ygr,Zdaughter,          &
                             Sdaughter_3D,rg_3D,i_HAx,j_HAy,   &
                             Rnearest_HA,Znearest_HA,dx,dy,    &
!*                             PI_nps,PJ_nps,DROT_nps,       &
!*                             YY,MO,DD,HH,MM,SS,               &
                             a,b,zref1,zref2,avehpbl,Sdaughter_3Dp)
!**********************************************************
!* author: Xingxiu Deng, Jan. 17, 2003
!*-objective:
!*  This SUBROUTINE is to calculate sharing factors for the 
!*  entire domain using the mother-daughter approach.
!*-Algorithm:
!*  Given a grid point, nine mothers may have their daughters 
!* (sharing factors) at that point, only the strongest one 
!* is to be kept at that point.
!*              . . . 
!*              . + . 
!*              . . . 
!*-For more information about the mother-daughter approach, please
!* refer the following paper:
!* Deng, X. and R. Stull, 2005: A mesoscale analysis method for surface
!* potential temperature in mountainous and coastal terrain. Mon. Wea.
!* Rev., 133, 389-408.
!**********************************************************
!*-Modification history:
!* - added ncycl to memorize how many pass for the point source
!*   travel to each grid point. (Jan. 23, 2003)
!* - added memorizing minimum travelling distance(r_min) for the point 
!*   source to each grid point. (Jan. 27, 2003)
!* - Considered real mc2 grid as grid points and obs as Honoured 
!*   Ancestor (first mother) (Feb. 20, 2003)
!* - added ibegin,iend,jbegin and jend to define the true range of 
!*   grid points for updating Sdaughter, in order to reduce CPU time 
!*   for all iterations (May. 2003)
!* - added refinement for coastal terrain (June. 2003)
!* - added influence region (R_inf) to limit ibegin,iend,jbegin and jend,
!*   again in order to reduce CPU time. The influence region should
!*   be larger than the influence region used in data analysis
!*   process (June, 2004)
!* - added refinement for mountain-top observations (July, 2004)
!* - removed lat/lon (HAlat/HAlon, latgr/longr) for source points
!*   and other grid points, remove PI_nps.... ( Jan. 23, 2006)
!* - increased the range of xgrmin,xgrmax,ygrmin,ygrmax for a possible
!*   larger domain (Jan. 23, 2004)
!* - added dy in order to work for different dx and dy (Jan. 24, 2006)
!* - added if_limit to check whether or not apply influence region when
!*   calculating sharing factors (0: no; 1: yes) (Jan. 24, 2006)
!* - removed the two refinements temporarily (Jan. 26, 2006)
!**********************************************************
         IMPLICIT none

!* input variables:
!* obs
         INTEGER :: n_HA
         CHARACTER (len=5) name_HA(n_HA)
         REAL :: HAx(n_HA),HAy(n_HA)
!*         REAL :: HAlat(n_HA),HAlon(n_HA)
         REAL :: z_HA(n_HA)
         REAL :: S_HA(n_HA)   ! always = 1.
!* grid
         INTEGER :: ni,nj
!*	 REAL :: latgr(ni,nj), longr(ni,nj)
         REAL :: xgr(ni,nj),ygr(ni,nj)
         REAL Zdaughter(ni,nj)

         REAL dx, dy
!*         REAL :: PI_nps, PJ_nps, DROT_nps
!*         INTEGER :: YY,MO,DD,HH,MM,SS

!* free parameters
         REAL :: a,b
         REAL :: zref1, zref2
         REAL :: avehpbl
             
!* Misc variables
         INTEGER :: ista, i,j, ii,jj,i1,j1,i2,j2,i0,j0
         INTEGER :: k,kk
         LOGICAL :: dbg
         PARAMETER (dbg=.false.)
         REAL :: xgrmin,xgrmax,ygrmin,ygrmax

!* intermediate variables:
!*************************************************
!* index_hax: x-dir index of Honoured ancestor's nearest grid point
!* index_hay: y-dir index of Honoured ancestor's nearest grid point
!* rmin: distance between the Honoured ancestor and its nearest grid point
!* S_mother: sharing factor for a temporary mother
!* Z_mother: terrain height for a temporary mother
!* ipt, jpt: index of the grid point to the south west of Honoured Ancestor
!* ireturn: indicator (=0, HA is inside the grid domain)
!*                    (<0, HA is outside the grid domain)
!**************************************************
         INTEGER :: index_hax, index_hay
         REAL :: rmin
         REAL :: S_mother
         REAL :: Z_mother
         INTEGER :: ipt, jpt
         INTEGER :: ireturn

!*************************************************
!* icycle: counter of iterative cycle
!* ibegin, iend, jbegin and jend: true range of grid points to update Sdaughter.
!*       they are depend on icycle. Added to reduce CPU time for all iterations. 
!* ibegin(icycle) = index_HAx - icycle
!* iend(icycle) = index_HAx + icycle
!* jbegin(icycle) = index_HAy - icycle
!* jend(icycle) = index_HAy + icycle
!* R_inf: Radius of influence region (m), used to 
!*        limit ibegin,iend,jbegin,jend. It should be larger than
!*        the radius of influence region used in data analysis.
!*        I.e., in ADAS if xyrange=90.km, cut-off radii in ADAS is 273.14km.
!*        R_inf = 300 km
!* if_limit: logical variable, to tell whether or not apply influence region
!*           when calculating sharing factor
!* tmp1,tmp2,tmp3,tmp4: temporary array
!* To save space, no dimension is used,
!*************************************************
         INTEGER :: icycle
         INTEGER :: ibegin, iend, jbegin, jend
         INTEGER :: ibegin_min, iend_max, jbegin_min, jend_max
         REAL :: tmp1,tmp2,tmp3,tmp4
         REAL R_inf
         PARAMETER (R_inf=300.E03)
         LOGICAL if_limit
         PARAMETER (if_limit=0)

!*************************************************
!* s: sharing factor from every possible mothers for a grid point
!* r: distance from every possible mothers to a grid point
!* imother: i index for every possible mothers 
!* jmother: j index for every possible mothers 
!*************************************************
         REAL :: s(9), r(9), imother(9),jmother(9)

!*************************************************
!* variables to memorize no. of pass and minimum traveling distance
!*************************************************
!* ncycl, ncycl0: variable to memorize no. of cycle for each grid point
!* r_min,rg: minimum distance from source point to the grid point
!*          (minimum of r_smax) (m)
!* rg0 :    temporary array
!* rg_out, rg_out0: similar to rg and rg0, but initialized to be 
!*          99999.E03(ini_rg) rather than 0.
!*          and final output rg_3D = rg_out rather than rg.
!* r_smax: distance from the source point to the grid point via 
!*         every possible mothers that has maximum contributions 
!*         (for one cycle only)
!*************************************************
         INTEGER ncycl(ni,nj),ncycl0(ni,nj)
         REAL rg(ni,nj),rg0(ni,nj)
         REAL rg_out(ni,nj),rg_out0(ni,nj)
         REAL r_min(ni,nj)
         REAL r_smax(ni,nj,9)
         REAL :: ini_rg       ! m
         PARAMETER (ini_rg = 99999.E03)     

!*************************************************
!* gx,gy:  x and y coord for the grid point considered (m)
!* R_gaussian: Gaussian e-folding distance=xyrange(1) in ADAS
!* Smother: sharing factor for every possible mother
!* Sdaughter: sharing factor for daughters
!* Stmp: temporary array for Sdaughter, used for comparison
!*       in iterative process.
!*************************************************
         REAL :: gx, gy
         REAL :: R_gaussian      ! m
         PARAMETER (R_gaussian = 90.E03)     
         REAL Smother(ni,nj)
         REAL Sdaughter(ni,nj)
         REAL Stmp(ni,nj)

!*************************************************
!* output variables 
!*
!* Sdaughter_3D(i,j,ista): sharing factor of grid point(i,j) 
!*           in response to HA at ista
!* Sdaughter_3Dp(i,j,ista): same as Sdaughter_3D but after Gaussian dropoff 
!* ncycl_3D: No. of pass from source - ista to every (i,j)
!*           = 0 if source does not travel to the point
!* rg_3D(i,j,ista): circuitous travel distance from HA at ista 
!*           to every grid point (i,j)
!*           = ini_rg(9999.E03m) if source does not travel to the point
!* following variables are output for use in ADAS obscor_new,brtobs3d_new
!* i_HAx: x-index of HAs' nearest grid point
!* j_HAy: y-index of HAs' nearest grid point
!* Rnearest_HA: distance between HA's nearest grid point and HA
!* Znearest_HA: height diff between HA's nearest grid point and HA
!*************************************************
         REAL Sdaughter_3D(ni,nj,n_HA),Sdaughter_3Dp(ni,nj,n_HA)
         INTEGER ncycl_3D(ni,nj,n_HA)
         REAL rg_3D(ni,nj,n_HA)
         INTEGER :: i_HAx(n_HA)
         INTEGER :: j_HAy(n_HA)
         REAL :: Rnearest_HA(n_HA)
         REAL :: Znearest_HA(n_HA)

!*---------- Begin executable code  ----------------
!* Initialize arrays
          Sdaughter_3D = 0.
          rg_3D = ini_rg
          ncycl_3D = 0
          Rnearest_HA = 0.
          Znearest_HA = 0.
          i_HAx = 0
          j_HAy = 0
          S_HA = 1.
!* find xgrmin,xgrmax,ygrmin,ygrmax
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

!* loop over all the Honoured Ancestors (observations)
         DO 1000 ista = 1, n_HA

!* Initialize the following two arrays for each HA
          Smother = 0.
          Sdaughter = 0.

!***********************************************
!* When Honoured Ancestor is not co-located with a grid point, 
!* an approximation is applied (the Honoured Ancestor is approximated
!* to be collocated with the grid point that is nearest
!* to the Honoured Ancestor in the sense of minimum elevation 
!* difference. The sharing factor to the nearest is reduced by 
!* elevation difference using a Gaussian drop-off with a standard
!* deviation of avehpbl.
!***********************************************

!* find grid point to the south and west of the HA first.
        IF ( HAx(ista) .LT. xgrmin .OR. HAx(ista) .GT. xgrmax .OR.    &
             HAx(ista) .LT. xgrmin .OR. HAx(ista) .GT. xgrmax  ) THEN
          ireturn = -1
!*         print*,'xgrmin,xgrmax,ygrmin,ygrmax',xgrmin,xgrmax,ygrmin,ygrmax
!*         print*,'HAx,HAy',HAx(ista),HAy(ista)
        ELSE
         CALL findlc(ni,nj,xgr,ygr,HAx(ista),HAy(ista),ipt,jpt,ireturn)
          print*, 'ipt,jpt,ireturn',ipt,jpt,ireturn
        ENDIF
        IF (ireturn .LT. 0 ) THEN
         print*,'Observation ', name_HA(ista), ' is outside the domain'
         GO TO 1000
        ENDIF
!* Find nearest point to the Honoured Ancestor if it is not collcated
!* with a grid point
        IF (HAx(ista) .NE. xgr(ipt,jpt) .AND.  &
            HAy(ista) .NE. ygr(ipt,jpt)) THEN

          CALL find_index_HA(ipt,jpt,HAx(ista),HAy(ista),Z_HA(ista),   &
                                  ni,nj,xgr,ygr,Zdaughter,             &
                                  index_HAx,index_HAy,rmin)
          i_HAx(ista) = index_HAx
          j_HAy(ista) = index_HAy
          Rnearest_HA(ista) = rmin
          Znearest_HA(ista) = Z_HA(ista)-Zdaughter(index_hax,index_hay)

!* Approximate HA to be collocated with the nearest neighboring grid point
!* with minimum elevation difference between them.
!* Need to account for elevation difference between HA and its nearest grid point.
!* avehpbl is approximated average PBL depth read from namelist &mdpara. If as follows, 
!* then Z_HA should be corrected to Zdaughter(index_hax,index_hay)) 
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
          print*, 'rmin,index_HAx,index_HAy',rmin,index_HAx,index_HAy
          print*, 'Z_HA', Z_HA(ista)
          Smother(index_HAx,index_HAy) = S_HA(ista)
          Sdaughter(index_HAx,index_HAy) = Smother(index_hax,index_hay)
!*          Z_HA(ista) = Z_HA(ista)
         ENDIF

!* Memorize sharing factors at the begining for comparison with those 
!* after one iteration
         Stmp = Sdaughter

!* 1st cycle 
         icycle = 1

!* Initialize ncycl and ncycl0 for every grid point
         ncycl = 0
         ncycl0 = ncycl
!* Initialize rg and rg0 for every grid point
         rg = 0.
         rg0 = rg
         rg_out = ini_rg
         rg_out(index_hax,index_hay) = 0.
         rg_out0 = rg_out

!* initialize ibegin, iend, jbegin, jend
          ibegin_min = index_HAx 
          iend_max = index_HAx 
          jbegin_min = index_HAy 
          jend_max = index_HAy 

          ibegin = ibegin_min
          iend = iend_max
          jbegin = jbegin_min
          jend = jend_max

10       CONTINUE
         IF (dbg) THEN
            print*,'icycle= ', icycle
         ENDIF
!* Given a grid point, nine mothers may have their daughters (sharing factors)
!* at that point, only the strongest one is to be kept at that point.
!* initialize sharing factors and distance from every possible mothers 
!* given a grid point for every cycle
         s = 0.
         r = 0.

!* Treat every grid point as a daughter, and find contributions from
!* its 8 neighbouring mothers and 1 mother of itself, and leave the 
!* largest one among 9 contributions.
!*         DO i = 1, ni
!*         DO j = 1, nj
!* Given HAx and HAy, we get index_HAx, index_HAy, we can update
!* Sdaughter in a gradually increasing areas rather than the entire domain
!* at the very begining. Use ibegin, iend, jbegin, jend instead.

!* Further define region of influence R_inf to limit
!* ibegin, iend, jbegin, jend
         IF (if_limit .eq. 1) THEN
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
          DO i = ibegin, iend
          DO j = jbegin, jend
         
          gx = xgr(i,j)
          gy = ygr(i,j)
! sharing factor if itself is mother
         s(9) = Sdaughter(i,j)
         r(9) = 0.
         imother(9) = i
         jmother(9) = j

! sharing factor if the point to the north-east is its mother
         i1 = i+1
         j1 = j+1
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(1) = i1
          jmother(1) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,1,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,1,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(1) = 0.
          r(1) = 0.
          imother(1) = i
          jmother(1) = j
         ENDIF

! sharing factor if the point to the east is its mother
         i1 = i+1
         j1 = j
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(2) = i1
          jmother(2) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,2,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,2,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(2) = 0.
          r(2) = 0.
          imother(2) = i
          jmother(2) = j
         ENDIF

! sharing factor if the point to the south-east is its mother
         i1 = i+1
         j1 = j-1
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(3) = i1
          jmother(3) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,3,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,3,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(3) = 0.
          r(3) = 0.
          imother(3) = i
          jmother(3) = j
         ENDIF
! sharing factor if the point ot the south is its mother
         i1 = i
         j1 = j-1
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(4) = i1
          jmother(4) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,4,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,4,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(4) = 0.
          r(4) = 0.
          imother(4) = i
          jmother(4) = j
         ENDIF
! sharing factor if the point to the south-west is its mother
         i1 = i-1
         j1 = j-1
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(5) = i1
          jmother(5) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,5,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,5,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(5) = 0.
          r(5) = 0.
          imother(5) = i
          jmother(5) = j
         ENDIF
! sharing factor if the point to the west is its mother
         i1 = i-1
         j1 = j
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(6) = i1
          jmother(6) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,6,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,6,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(6) = 0.
          r(6) = 0.
          imother(6) = i
          jmother(6) = j
         ENDIF
! sharing factor if the point to the north-west is its mother
         i1 = i-1
         j1 = j+1
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(7) = i1
          jmother(7) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,7,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,7,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(7) = 0.
          r(7) = 0.
          imother(7) = i
          jmother(7) = j
         ENDIF
! sharing factor if the point to the north is its mother
         i1 = i
         j1 = j+1
         IF ( i1 .GE. 1 .and. i1 .LE. ni            &
           .and. j1 .GE. 1 .and. j1 .LE. nj )  THEN
          imother(8) = i1
          jmother(8) = j1

          CALL mother(ista,z_HA,n_HA,ni,nj,Smother,i1,j1,s,8,Zdaughter,i,j, &
                      a,b,zref1,zref2)
          CALL distance(ista,i1,j1,r,8,i,j,gx,gy,xgr,ygr,ni,nj)
         ELSE
          s(8) = 0.
          r(8) = 0.
          imother(8) = i
          jmother(8) = j
         ENDIF
! Keep the strongest daughter at the given point (i,j)        
         Sdaughter(i,j) = amax1(s(9),s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8))

! Find the nearest mother that has maximum contribution to the grid point
! r_smax: the distance from source point to the grid point via the 
! the mother point that has largest contribution.
! rg0(i2,j2): the distance from source point to the nearest mother point
! r(k): distance from grid point to the mother point that has largest 
! contribution. 
!                     rg0(i2,j2)            r(k)
!  HA(nearest point) ------------> mother -------> grid point
!
! Sometimes, there are several possible mothers that have same contribution
! to a grid point, the minimum distance is considered. 
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
! ELSE
! no changes to travelling dictance
      ENDIF

      ENDDO    ! end of i
      ENDDO    ! end of j

! Till now, one cycle is finished for every grid point (i,j)
! Following lines are for debug and to get ready for next cycle

! To memorize the maximum no. of cycle for every gridpoint
! and memorize the distance from the Honoured Ancestor to 
! every grid point 
         Do i =1, ni
         DO j = 1, nj
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

!  Update Smother after one cycle
         Smother = Sdaughter

!  update ncycl0, and rg0
         ncycl0 = ncycl
         rg0 = rg
         rg_out0 = rg_out
        
! Check if Sdaughter over entire domain update against to
! the values of the previous cycle

         DO i =1, ni
         DO j =1, nj
          IF (abs(Sdaughter(i,j)-Stmp(i,j)) .GE. 0.01) GO TO 888
         ENDDO
         ENDDO
! No change between current and previous cycle for all points, exit cycle
         GO TO 999

! Yes, find change between current and previous cycle
! memorize current values and go to next cycle
888      CONTINUE

         IF (dbg) THEN
           WRITE(60,*) 'icycle', icycle
           DO j0=nj,1,-1
            WRITE(60,110) (Sdaughter(i0,j0),i0=1,ni)
           ENDDO
           print*,'ncycl'
           DO j0=nj,1,-1
            WRITE(6,111) (ncycl(i0,j0),i0=1,ni)
           ENDDO
 
           print*,'circuitous travel distance'
           DO j0=nj,1,-1
            WRITE(6,112) (rg(i0,j0)/1000.,i0=1,ni)    ! print in km
           ENDDO
         ENDIF
110      FORMAT(13f6.2)
111      FORMAT (13I5)
112      FORMAT(13f8.2)

         Stmp = Sdaughter
         icycle = icycle + 1 
         GO TO 10

! exit cycle
999      CONTINUE

       IF (dbg) THEN
        WRITE(99,*) ista,index_HAx,index_HAy
        WRITE(99,*) ibegin, iend, jbegin, jend 
        WRITE(99,*) ibegin_min, iend_max, jbegin_min, jend_max 
        WRITE(99,*) tmp1,tmp2,tmp3,tmp4
        WRITE(99,*) 'cycle= ', icycle
       ENDIF

!* for plotting purpose (plot sharing factor after Gaussian drop-off)
!* should comment following for analysis purpose,as this will be 
!* done in ADAS Bratseth scheme. Now use a new variable (Mar. 3, 2006)
!* Adding consideration of horizontal distance using Gaussian function:
!* rg_out (meters), Rnearest_HA(meters), R_gaussian(meters)
         DO j0 = 1, nj
          DO i0 = 1, ni
           Sdaughter_3Dp(i0,j0,ista) = Sdaughter(i0,j0) *     &
              exp(-0.5*((rg_out(i0,j0)+Rnearest_HA(ista))/R_gaussian)**2)
           ENDDO
         ENDDO

         Do j =1,nj
         Do i =1,ni
          IF (abs(Zdaughter(i,j)-Z_HA(ista)).GE.Zref2) rg_out(i,j)= ini_rg ! m
         ENDDO
         ENDDO

         Do j =1,nj
         Do i =1,ni
          Sdaughter_3D(i,j,ista) = Sdaughter(i,j)
          ncycl_3D(i,j,ista) = ncycl(i,j)
          rg_3D(i,j,ista) = rg_out(i,j) + Rnearest_HA(ista)      ! meters
         ENDDO
         ENDDO

 1000    CONTINUE       ! ista

!!*        CALL output_3D(Sdaughter_3Dp,Zdaughter,ncycl_3D,rg_3D/1000., &
!!*                     ni,nj,n_HA,  &
!!*                     name_HA,latgr,longr,  &
!!*                     PI_nps,PJ_nps,DROT_nps,dx,YY,MO,DD,HH,MM,SS)
!!*
!!*        CALL output_HA(n_HA,name_HA,HAlat,HAlon,YY,MO,DD,HH,MM,SS) 
!!* to plot sharing factor for only one station and for idealized domain
!*         CALL output(Sdaughter_3Dp, Zdaughter,ncycl_3D,  &
!*                     rg_3D/1000.,ni,nj,n_HA)
!*         CALL output_obs(n_HA,name_HA,HAx,HAy,z_HA,S_HA)

         RETURN
         END

         SUBROUTINE mother(ista,z_HA,n_HA,ni,nj,Smother,im,jm,ss,in,   &
                             Zdaughter,id,jd,a,b,zref1,zref2)
!***********************************************************
!* author: Xingxiu Deng, Jan. 20, 2003
!* This SUBROUTINE is to clculate sharing factors of a grid(id,jd) 
!* point from one of its mother (Smother(im,jm))
!***********************************************************
         IMPLICIT NONE
!*-----------------------------------------------------------
!* input variables
!* ni : x-dir dimension size
!* nj : y-dir dimension size
!* ista: index of the Honoured Ancestor
!* im,jm: grid index of a mother
!* id,jd: grid index of a daughter
!* in: order (or sequence) of the mother around the grid point
!* Z_mother : elevation of the mother (m)
!* Zdaughter: elevation of the daughter (m)
!* Smother : sharing factors of mothers
!*-----------------------------------------------------------
!* output variable
!* ss :  sharing factors at the point(id,jd) from the mother (im,jm)
         INTEGER :: n_HA
         REAL :: z_HA(n_HA)
         REAL :: a,b,zref1,zref2

         INTEGER :: ni,nj,ista,im,jm,id,jd,in
         REAL :: Z_mother
         REAL :: Zdaughter(ni,nj),Smother(ni,nj)
         REAL :: ss(9)
      
! misc variables
         INTEGER :: i,j

         Z_mother = Zdaughter(im,jm)

         ss(in) = Smother(im,jm) *                        &
          (1.-((abs(Z_mother - Zdaughter(id,jd)))/zref1)**a) *  &
          (1.-((abs(z_HA(ista) - Zdaughter(id,jd)))/zref2)**b)

         IF ( abs(z_HA(ista) - Zdaughter(id,jd)) .GT. zref2 ) then
           ss(in) = 0.0
         ENDIF
         IF ( abs(Z_mother - Zdaughter(id,jd)) .GT. zref1 ) then
           ss(in) = 0.0
         ENDIF
 
         RETURN
         END

         SUBROUTINE distance(ista,im,jm,r,in,    &
                             id,jd,gx,gy,xgr,ygr,ni,nj)
!**********************************************************
!* author: Xingxiu Deng, Jan. 27, 2003
!*   This SUBROUTINE is to calculate distance of a grid(id,jd) 
!*   point to one of its mother (im,jm)
!**********************************************************
         IMPLICIT NONE
         
!*-----------------------------------------------------------
! input variables
!* ista: index of the Honoured Ancestor
!* im,jm: grid index of a mother
!* id,jd: grid index of a daughter
!* gx,gy: x and y coord of a grid point (m)
!* in: order (or sequence) of the mother around the grid point
!*-----------------------------------------------------------
!* output variable
!* r : distance from the point(id,jd) to the mother (im,jm)
!**********************************************************
         INTEGER :: ista,im,jm,id,jd,in
         REAL :: r(9)
      
         REAL :: gx,gy    ! x and y coord for the grid point (m)
         REAL :: mx,my    ! x and y coord for the mother (m)
         INTEGER :: ni,nj ! grid dimension size
         REAL :: xgr(ni,nj),ygr(ni,nj) ! x and y coord of grid points
         
         mx = xgr(im,jm)
         my = ygr(im,jm)
         
         r(in) = SQRT((gx -mx)**2 + (gy-my)**2) 

         RETURN
         END

        SUBROUTINE find_index_HA(ipt,jpt,xpt,ypt,zpt,nx,ny,xg,yg,zg,   &
                                  index_xpt,index_ypt,rmin)
!**********************************************************
!* author: Xingxiu Deng, Feb. 20, 2003
!*   This SUBROUTINE is to find a nearest grid point to the obs point
!*   (xpt,ypt) given the grid point to the south and west (ipt,jpt)
!* Modifing history
!*   March 05, 2003: changed to find the nearest grid point based on 
!*                   elevation difference rather than horizontal 
!*                   distance. rmin now is defined as the horizontal
!*                   distance between HA and the nearest grid point.
!**********************************************************
        IMPLICIT NONE
         
        INTEGER :: ipt     ! i index to the west of desired
                           ! location (xpt, ypt)
        INTEGER :: jpt     ! j index to the south of desired
                           ! location (xpt,ypt)
        REAL :: xpt, ypt,zpt ! x,y,z coord of desired location
         
        INTEGER :: nx,ny   ! Grid dimensions 
        REAL :: xg(nx,ny)  ! x coordinate of grid points in
                           ! physical/comp. space (m)
        REAL :: yg(nx,ny)  ! y coordinate of grid points in
                           ! physical/comp. space (m)
        REAL :: zg(nx,ny)  ! z coordinate of grid points in 
                           ! physical space (m)

        INTEGER :: index_xpt   ! i index that is closest to the desired
                               ! location (xpt, ypt)
        INTEGER :: index_ypt   ! j index that is closest to the desired
                               ! location (xpt, ypt)

        REAL :: r(4),zr(4)           ! intermediate array
        REAL :: rmin, zrmin
        INTEGER :: k, kk,kz

!*---------- Begin executable code ----------------------
         r = 0.

         r(1) = SQRT ( (xg(ipt,jpt) - xpt)**2 +        &
                       (yg(ipt,jpt) - ypt)**2 )
         zr(1) = SQRT ( (zg(ipt,jpt) - zpt)**2 )
         r(2) = SQRT ( (xg(ipt,jpt+1) - xpt)**2 +      &
                       (yg(ipt,jpt+1) - ypt)**2 )
         zr(2) = SQRT ( (zg(ipt,jpt+1) - zpt)**2)
         r(3) = SQRT ( (xg(ipt+1,jpt+1) - xpt)**2 +    &
                       (yg(ipt+1,jpt+1) - ypt)**2 ) 
         zr(3) = SQRT ( (zg(ipt+1,jpt+1) - zpt)**2 )
         r(4) = SQRT ( (xg(ipt+1,jpt) - xpt)**2 +      &
                       (yg(ipt+1,jpt) - ypt)**2  )
         zr(4) = SQRT ( (zg(ipt+1,jpt) - zpt)**2)
         
!         rmin = amin1 ( r(1), r(2), r(3),r(4) )
!         DO k = 1, 4
!           IF ( rmin .eq. r(k) ) kk = k
!         ENDDO
!         IF ( kk .eq. 1 ) THEN
!          index_xpt = ipt
!          index_ypt = jpt
!         ELSE IF ( kk .eq. 2 ) THEN
!          index_xpt = ipt
!          index_ypt = jpt + 1
!         ELSE IF ( kk .eq. 3 ) THEN
!          index_xpt = ipt + 1
!          index_ypt = jpt + 1
!         ELSE IF ( kk .eq. 4 ) THEN
!          index_xpt = ipt + 1
!          index_ypt = jpt 
!         ENDIF

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
         rmin = r(kz)
         print*, 'rmin,index_xpt,index_ypt',rmin,index_xpt,index_ypt
         print*, 'zg', zg(index_xpt,index_ypt)
         RETURN
         END

!************************************
        SUBROUTINE output_HA(n_HA,name_HA,HAlat,HAlon,YY,MONTH,DD,HH,MM,SS) 
        IMPLICIT none

        INTEGER n_HA
        CHARACTER(len=5):: name_HA(n_HA)
        REAL :: S_HA(n_HA)               ! sharing factor for each station
        REAL :: HAlat(n_HA),HAlon(n_HA)  ! lat/lon of obs station
        INTEGER :: YY,MONTH,DD,HH,MM,SS
        character (len=3) :: mon, mo         ! month
        character (len=2) :: cday            ! day
        character (len=4) :: cyear           ! year

        CHARACTER(len=40):: stringZ

       INTEGER maxHA
       PARAMETER (maxHA = 500)
       CHARACTER(len=8):: stid(maxHA)

       REAL :: tim
       INTEGER :: nlev
       INTEGER :: nflag

       INTEGER :: ista,ivar,ik, i

! as S_HA is always 1.0, we set its values here rather through virtue variable
       S_HA = 1.
! convert integer day to  character cday
      WRITE(cday,'(I2)') DD
      WRITE(cyear,'(I4)') YY
      print*, 'cday',cday
      print*, 'cyear',cyear
      print*, mo(MONTH)

! for DD < 10
      DO i=1, 2
       IF ( cday (i:i) .eq. ' ') cday(i:i)='0'
      ENDDO

       tim = 0.0
       nlev = 1
       nflag = 1

! Get stid with length of 8 from nme_HA of length of 5
! stid must have a length of 8 in order for it to be a station file
       DO ista = 1, n_HA
       stid(ista) = name_HA(ista)//'   '
       ENDDO

! write out Grads data
        OPEN (61,file='HA_grads.dat',status='unknown',form='unformatted',  &
              access='sequential')
        DO ista=1, n_HA
!sequence: stid(character(len=8)),lat,lon,tim,nlev,nflag
        WRITE(61) stid(ista),HAlat(ista),HAlon(ista),tim,nlev,nflag
        WRITE(61) S_HA(ista),-999.
        ENDDO
         ivar = 2
         ik = 0
! Tell GrADS the end of record
        NLEV = 0
        WRITE(61) stid(n_HA),HAlat(n_HA),HAlon(n_HA),tim,nlev,nflag
        CLOSE(61)

! open GrADS descriptor file
       OPEN (62, file='HA_grads.ctl', form='formatted',  &
                 STATUS = 'unknown')
       WRITE(62,'(''DSET'',a30)') 'HA_grads.dat'
       write(62,'(''byteswapped'')')
       WRITE(62,'(''DTYPE'',a30)') 'station'
       WRITE(62,'(''STNMAP'',a30)') 'HA_grads.map'
       write(62,'(''OPTIONS sequential'')')
       write(62,'(''UNDEF'',F12.6)') -999.
       write(62,108) 1, HH, cday//mo(month)//cyear, 1, 'HR'
       write(62,109)  ivar
       write(62,110) 'S_HA  ',ik,99,'Sharing factor of the Honoured Ancestor'
       write(62,110) 'dummy ',ik,99,'dummy variable with -999.'
       write(62,111) 'ENDVARS'
108    format('TDEF ',I2,' LINEAR ',I2,'Z',A9,1X,I2,A2)
109    format('VARS ',I3)
110    format(A6,I3,I3,1X,A40)
111    format(A7)
       CLOSE(62)

        RETURN
        END

      SUBROUTINE output_3D(Sd,Zd,ncycl,rg,ni,nj,n_HA,name_HA,latgr,longr, &
                          PI_nps,PJ_nps,DROT_nps,dx,YY,MONTH,DD,HH,MM,SS)
        IMPLICIT none
        INTEGER :: ni,nj,n_HA
        CHARACTER (len=5) name_HA(n_HA)
        INTEGER :: int_name(n_HA)
        REAL :: Sd(ni,nj,n_HA), Zd(ni,nj)
        REAL :: rg(ni,nj,n_HA)
        INTEGER :: ncycl(ni,nj,n_HA)
 
        REAL :: latgr(ni,nj), longr(ni,nj)
        REAL :: PI_nps,PJ_nps,DROT_nps,dx
        REAL :: lamin,lamax,lomin,lomax,deltala,deltalo
        INTEGER :: YY,MONTH,DD,HH,MM,SS
        character (len=3) :: mon, mo         ! month
        character (len=2) :: cday         ! day
        character (len=4) :: cyear        ! year

        CHARACTER(len=40):: stringZ
        INTEGER :: ivar, irec

        INTEGER :: i,j,k

! get int_name from name_HA
      DO k = 1, n_HA
       read(name_HA(k)(2:5),'(i4)') int_name(k)
      ENDDO

! find deltala, deltalo for PDEF
      lamin = latgr(1,1)
      lamax = latgr(1,1)
      lomin = longr(1,1)
      lomax = longr(1,1)
      DO j =1, nj
      DO i =1, ni
       if ( lamin .gt. latgr(i,j) ) lamin = latgr(i,j)
       if ( lamax .lt. latgr(i,j) ) lamax = latgr(i,j)
      ENDDO
      ENDDO

      DO j =1, nj
      DO i =1, ni
       if ( lomin .gt. longr(i,j) ) lomin = longr(i,j)
       if ( lomax .lt. longr(i,j) ) lomax = longr(i,j)
      ENDDO
      ENDDO
      deltalo = (lomax - lomin)/(ni-1)
      deltala = (lamax - lamin)/(nj-1)

      print *, 'Deng check lamax',lamax,lamin,lomax,lomin
      print *, 'Deng check PI_nps', PI_nps, PJ_nps, DROT_nps, dx

!   convert integer day to  character cday
      WRITE(cday,'(I2)') DD
      WRITE(cyear,'(I4)') YY
      print*, 'cday',cday
      print*, 'cyear',cyear
      print*, mo(MONTH)

! for DD < 10
      DO i=1, 2
       IF ( cday (i:i) .eq. ' ') cday(i:i)='0'
      ENDDO

! write out Grads data
      OPEN (60,file='Sdaughter.dat',status='unknown',form='unformatted',  &
              access='direct', recl = 4*ni*nj)
         irec = 1
         ivar = 1
         DO k=1,n_HA
         write(60,rec=irec) ((Sd(i,j,k), i=1,ni),j=1,nj)
         irec = irec + 1
         ENDDO
         ivar = ivar + 1
         write(60,rec=irec) ((Zd(i,j), i=1,ni),j=1,nj)
         irec = irec + 1
         ivar = ivar + 1
         DO k=1,n_HA
         write(60,rec=irec) ((float(ncycl(i,j,k)), i=1,ni),j=1,nj)
         irec = irec + 1
         ENDDO
         ivar = ivar + 1
         DO k=1,n_HA
         write(60,rec=irec) ((rg(i,j,k), i=1,ni),j=1,nj)
         irec = irec + 1
         ENDDO
         CLOSE(60)
         
! write ctl file 
         open(70, file='Sdaughter.ctl', status='unknown', form='formatted')
         write(70,'(''DSET'',a30)') 'Sdaughter.dat'
         write(70,'(''OPTIONS little_endian'')')
         write(70,'(''byteswapped'')')
         write(70,'(''UNDEF'',F12.6)') -999.
         write(70,104) ni,lomin,deltalo
         write(70,105) nj,lamin,deltala
         write(stringZ,'(a,i5,a)') '(''ZDEF'', I4, '' LEVELS''', n_HA,'f10.2)'
         write(70,stringZ) n_HA,(float(int_name(k)),k=1,n_HA) !(float(k),k=1,n_HA)
         write(70,107) ni,nj, PI_nps,PJ_nps,DROT_nps, dx/1000.
         write(70,108) 1, HH, cday//mo(month)//cyear, 1, 'HR'
         write(70,109)  ivar      ! No. of variables
         write(70,110) 'Sd    ', n_HA, 99, ' sharing factors'
         write(70,110) 'Zd    ', 1, 99, ' terrain height'
         write(70,110) 'ncycl ', n_HA, 99, ' No. of pass   '
         write(70,110) 'rmin  ', n_HA, 99, ' traveling distance'
         write(70,111) 'ENDVARS'
104     format('XDEF ',I5,' LINEAR ',2F)
105     format('YDEF ',I5,' LINEAR ',2F)
107     format('PDEF ',I3,1X,I3,' nps ',4f12.6)
108     format('TDEF ',I2,' LINEAR ',I2,'Z',A9,1X,I2,A2)
109     format('VARS ',I3)
110     format(A6,I3,I3,1X,A40)
111     format(A7)
        CLOSE(70)
        RETURN
        END

      SUBROUTINE output(Sd, Zd,ncycl, rg,ni,nj,n_HA)
      IMPLICIT none
      INTEGER :: ni,nj,n_HA
      REAL :: Sd(ni,nj,n_HA), Zd(ni,nj)
      REAL :: rg(ni,nj,n_HA)
      INTEGER :: ncycl(ni,nj,n_HA)

      CHARACTER(len=40):: stringZ
      INTEGER :: ivar, irec

      INTEGER :: i,j,k

! write out Grads data
      OPEN (60,file='Sdaughter.dat',status='unknown',form='unformatted',  &
               access='direct', recl = 4*ni*nj)
      irec = 1
      ivar = 1
      DO k=1,n_HA
      write(60,rec=irec) ((Sd(i,j,k), i=1,ni),j=1,nj)
      irec = irec + 1
      ENDDO
      ivar = ivar + 1
      write(60,rec=irec) ((Zd(i,j), i=1,ni),j=1,nj)
      irec = irec + 1
      ivar = ivar + 1
      DO k=1,n_HA
      write(60,rec=irec) ((float(ncycl(i,j,k)), i=1,ni),j=1,nj)
      irec = irec + 1
      ENDDO
      ivar = ivar + 1
      DO k=1,n_HA
      write(60,rec=irec) ((rg(i,j,k), i=1,ni),j=1,nj)
      irec = irec + 1
      ENDDO
      CLOSE(60)
         
! write ctl file 
      open(70, file='Sdaughter.ctl', status='unknown', form='formatted')
      write(70,'(''DSET'',a30)') 'Sdaughter.dat'
      write(70,'(''OPTIONS little_endian'')')
      write(70,'(''UNDEF'',F12.6)') -999.
      write(70,104) ni,1.0,1.0
      write(70,105) nj, 1.0, 1.0
      write(stringZ,'(a,i5,a)') '(''ZDEF'',I4,'' LEVELS''',n_HA,'f10.2)'
!*      write(70,stringZ) 1, 0.
      write(70,stringZ) n_HA,(float(k),k=1,n_HA)
      write(70,108) 1, 1, '20JAN2003', 1, 'HR'
      write(70,109)  ivar      ! No. of variables
      write(70,110) 'Sd    ', n_HA, 99, ' sharing factors'
      write(70,110) 'Zd    ', 1, 99, ' terrain height'
      write(70,110) 'ncycl ', n_HA, 99, ' No. of pass   '
      write(70,110) 'rmin  ', n_HA, 99, ' traveling distance'
      write(70,111) 'ENDVARS'
104   format('XDEF ',I5,' LINEAR ',2F)
105   format('YDEF ',I5,' LINEAR ',2F)
!107  format('PDEF ',I3,1X,I3,' nps ',4f8.2)
108   format('TDEF ',I2,' LINEAR ',I2,'Z',A9,1X,I2,A2)
109   format('VARS ',I3)
110   format(A6,I3,I3,1X,A40)
111   format(A7)
      CLOSE(70)
      RETURN
      END

!************************************
      function mo(mm) result (mon)
      character (len=3) :: mon
      INTEGER :: mm
      if(mm.eq.1)mon='JAN'
      if(mm.eq.2)mon='FEB'
      if(mm.eq.3)mon='MAR'
      if(mm.eq.4)mon='APR'
      if(mm.eq.5)mon='MAY'
      if(mm.eq.6)mon='JUN'
      if(mm.eq.7)mon='JUL'
      if(mm.eq.8)mon='AUG'
      if(mm.eq.9)mon='SEP'
      if(mm.eq.10)mon='OCT'
      if(mm.eq.11)mon='NOV'
      if(mm.eq.12)mon='DEC'
      end function mo
!********************************
       SUBROUTINE output_obs(n_HA,name_HA,HAx,HAy,z_HA,S_HA)
       IMPLICIT none

       CHARACTER(len=40):: stringZ

       INTEGER n_HA
       REAL HAx(n_HA),HAy(n_HA),z_HA(n_HA),S_HA(n_HA)
       CHARACTER(len=5):: name_HA(n_HA)
       CHARACTER(len=8):: stid(n_HA)

       REAL :: tim
       INTEGER :: nlev
       INTEGER :: nflag

       INTEGER :: ista,ivar,ik, hour

       tim = 0.0
       nlev = 1
       nflag = 1

! Get stid with length of 8 from nme_HA of length of 5
! stid must have a length of 8 in order for it to be a station file
       DO ista = 1, n_HA
       stid(ista) = name_HA(ista)//'   '
       ENDDO

! write out Grads data
       OPEN (61,file='Obs_grads.dat',status='unknown',form='unformatted', &
              access='sequential')
        DO ista=1, n_HA
!* sequence: stid(character(len=8)),lat,lon,tim,nlev,nflag
        WRITE(61) stid(ista),HAy(ista),HAx(ista),tim,nlev,nflag
        WRITE(61) S_HA(ista)
        ENDDO
         ivar = 1
         ik = 0
         hour = 1
! Tell GrADS the end of record
        NLEV = 0
        WRITE(61) stid(n_HA),HAy(n_HA),HAx(n_HA),tim,nlev,nflag
        CLOSE(61)

! open GrADS descriptor file
       OPEN (62, NAME='Obs_grads.ctl', form='formatted',  &
                 STATUS = 'unknown')
       WRITE(62,'(''DSET'',a30)') 'Obs_grads.dat'
       WRITE(62,'(''DTYPE'',a30)') 'station'
       WRITE(62,'(''STNMAP'',a30)') 'Obs_grads.map'
       write(62,'(''OPTIONS sequential'')')
       write(62,'(''UNDEF'',F12.6)') -999.
       write(62,108) 1, hour, '20JAN2003', 1, 'HR'
       write(62,109)  ivar
       write(62,110) 'S_HA  ',ik,99,'Sharing factor of the source points'
       write(62,111) 'ENDVARS'
108     format('TDEF ',I2,' LINEAR ',I2,'Z',A9,1X,I2,A2)
109     format('VARS ',I3)
110     format(A6,I3,I3,1X,A40)
111     format(A7)
        CLOSE(62)

        RETURN
        END
