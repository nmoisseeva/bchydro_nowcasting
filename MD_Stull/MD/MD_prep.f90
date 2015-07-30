      subroutine MD_prep(outfstpf,nobs,stname,latstn_rot,lonstn_rot,elevstn, &
                    xyzstn_rot_8,iswc,jswc,ireturnc,nxc,nyc,topoc,nx,ny,topo, &
                    latrot,lonrot,xyzgr_rot_8,Rxy,indir)
!*********************************************************************
! author:	Xingxiu Deng
! date:		March - April 2007
! Object:	to prepare for calling the main program MD
!*********************************************************************
! input:        
!       for surface stations:
!               nobs: No. of surface stations
!               stname: name or character id of stations
!               latstn_rot: rotated latitude of stations
!               lonstn_rot: rotated longitude of stations
!               elevstn: elevation of stations (m)
!               xyzstn_rot_8: rotated Cartesian coords(3) of station
!               iswc: i index of the GP to the south and west of stn (wrt coarse grid)
!               jswc: j index of the GP to the south and west of stn (wrt coarse grid)
!               ireturnc: < 0, the station is outside the coarse grid 
!                         iswc=-999 and jswc=-999.
!       for coarse grid: 
!               nxc,nyc: dimension of the coarse grid
!               topoc: topography of the coarse grid (m)
!       for fine grid (work grid for MD): 
!               nx,ny: dimension of the fine grid
!               topo: topography of the fine grid (m)
!               latrot: rotated lat
!               lonrot: rotated lon
!               xyzgr_rot_8: rotated Cartesian coords(3)
!
!               Rxy: Radius of influence (km)
!               indir: directory where holds inputs
! output:       outfstpf: path/file name of output sharing factor file in
!                         RPN standard format
! derived output:
!*********************************************************************
      implicit none
!*** input 
!* - for surface stations
      integer nobs
      character (LEN=8) :: stname(nobs)
      real latstn_rot(nobs),lonstn_rot(nobs)
      real :: elevstn(nobs)
      real*8 xyzstn_rot_8(nobs,3) ! rotated Cartesian coords(3) of station
      integer iswc(nobs),jswc(nobs)
      integer ireturnc(nobs)


!* - for coarse domain
      integer nxc, nyc
      real topoc(nxc,nyc)

!* - for fine domain (or working grid for MD)
      integer nx,ny
      real topo(nx,ny)
      real*8 xyzgr_rot_8(nx,ny,3)
      real latrot(ny), lonrot(nx)

      integer gfactor
!*      
      real Rxy
      character(len=70) indir

!*** output 
      character(len=100) outfstpf

      real,allocatable :: Sdaughter_obstn(:,:)  ! sharing factors from source point to sfc station
      real,allocatable :: rg_obstn(:,:)         ! travel distance from source point to sfc station
      real,allocatable :: corr_obstn(:,:)       ! correlation between source point and obstn 
                                                !  [G(s)*S=G(rg_obstn)*Sdaughter_obstn]
      real,allocatable :: topofmc(:,:)         ! topography difference (fine-coarse) at source points 

!*   intermediate output
      integer,allocatable :: iif_obstn (:)       ! i-index of the obstn's nearest GP within fine grid
      integer,allocatable :: jjf_obstn (:)       ! j-index of the obstn's nearest GP within fine grid
      real,allocatable :: Rnearest_obstn (:)    ! distance between the GP at (iif_obstn, jjf_obstn) and obstn
      real, allocatable :: Znearest_obstn (:)   ! height diff between the GP at (iif_obstn, jjf_obstn) and obstn

!*   intermediate output
      real,allocatable :: Sdaughter_3D(:,:,:)   ! sharing factors at working grid points from all source points
      real,allocatable :: rg_3D(:,:,:)          ! circuitous travel distance from source point to the GP
      integer,allocatable :: i_HAx (:)          ! i-index of the source point's nearest GP
      integer,allocatable :: j_HAy (:)          ! j-index of the source point's nearest GP
      real,allocatable :: Rnearest_HA (:)       ! distance between the GP at (i_HAx, j_HAy) and the source point
      real, allocatable :: Znearest_HA (:)      !  height diff between the GP at (i_HAx, j_HAy) and the source point
      real :: mothera, motherb, zref1, zref2,avehpbl  ! MD free parameters
      real,allocatable :: Sdaughter_3Dp(:,:,:)  ! same as Sdaughter_3D but after Gaussian dropoff
      character (LEN=40) runname

!*intermediate
      real :: dx, dy
      integer nx_sp,ny_sp,nxy_sp
      integer nxy_sp_max
      parameter (nxy_sp_max=4)
      real,allocatable :: topo_sp(:,:)
      real,allocatable :: latrot_sp(:), lonrot_sp(:)
      real*8,allocatable :: xyzsp_rot_8(:,:,:)

      integer nobsng
      character(LEN=5),allocatable :: stnsng(:)
      real*8,allocatable :: xyzsng_rot_8(:,:)
      real,allocatable :: latrotsng(:),lonrotsng(:)
      real,allocatable :: hgtsng(:)
      real,allocatable :: qualsng(:)

!*misc
!     iswf, jswf: i,j index wrt the fine grid of GP corresponding to iswc
!     IIC, JJC: i,j index wrt the coarse grid of station's neighbouring GPs
!     IIf, JJf:  i,j index wrt the fine grid corresponding to IIC and JJC
      integer ista,ijsp,i,j,k,iobsng,istatus
      integer IIc,JJc, IIf, JJf,II,JJ
      integer iswf,jswf
      logical :: dbg=.true.

      print*, '---- EXECUTE MD_prep -------------'

! --------------------------------------------------------------------
!*** calculate dx and dy for the working grid
      call cal_dist_sphere(xyzgr_rot_8(nx/2,ny/2,:),         &
                           xyzgr_rot_8(nx/2+1,ny/2,:),dx)
      print*, 'dx = ', dx
      call cal_dist_sphere(xyzgr_rot_8(nx/2,ny/2,:),         &
                           xyzgr_rot_8(nx/2,ny/2+1,:),dy)
      print*, 'dy = ', dy

!*** read  grid spacing ratio of the coarse to fine grid
      open(30,file=trim(indir)//'/gc2f_ratio',form='formatted',   &
              status='old')
      read(30,'(i)') gfactor 
      close(30)
      print*, 'gfactor= ', gfactor

!*** allocate output variables and intermediate output variables
        allocate(Sdaughter_obstn(nobs,nxy_sp_max),STAT=istatus)   
        allocate(rg_obstn(nobs,nxy_sp_max),STAT=istatus)        
        allocate(corr_obstn(nobs,nxy_sp_max),STAT=istatus)   
        allocate(topofmc(nobs,nxy_sp_max),STAT=istatus)

        allocate(iif_obstn(nobs),STAT=istatus)
        allocate(jjf_obstn(nobs),STAT=istatus)
        allocate(Rnearest_obstn(nobs),STAT=istatus)
        allocate(Znearest_obstn(nobs),STAT=istatus)

! --------------------------------------------------------------------
!*** For each surface station, treat its four neighbouring coarse grid points
!*   as source points and calculate sharing factors from each source point

      DO 2000 ista=1,nobs

       IF (ireturnc(ista) .LT. 0) THEN
        Sdaughter_obstn=-999.
        rg_obstn=-999.
        corr_obstn=-999.
        topofmc=-999.
        GO TO 2000
       ENDIF

        iswf=(iswc(ista)-1)*gfactor+1
        jswf=(jswc(ista)-1)*gfactor+1

        if (dbg) then
           print*, 'ista= ', ista
           print*, 'iswc, jswc= ', iswc(ista),jswc(ista)
           print*, 'iswf, jswf= ', iswf,jswf
        endif

        nx_sp=2
        ny_sp=2
        nxy_sp=nx_sp*ny_sp

        IF (nxy_sp .NE. nxy_sp_max) THEN
           print*, 'nxy_sp (',nxy_sp,') is not equal to nxy_sp_max (' &
                   ,nxy_sp_max,')'
           print*, 'redefine nxp_sp_max in MD_prep'
           stop
        ENDIF

        allocate(topo_sp(nx_sp,ny_sp),STAT=istatus)
        allocate(latrot_sp(ny_sp),lonrot_sp(nx_sp))
        allocate(xyzsp_rot_8(nx_sp,ny_sp,3),STAT=istatus)

        ijsp=0
        DO j=1,ny_sp
        DO i=1,nx_sp

        ijsp=ijsp+1

           IIf=(i-1)*gfactor+iswf
           JJf=(j-1)*gfactor+jswf
           IIc=iswc(ista)+i-1
           JJc=jswc(ista)+j-1

           DO k=1,3
              xyzsp_rot_8(i,j,k)=xyzgr_rot_8(IIf,JJf,k)
           ENDDO
           topo_sp(i,j)=topo(IIf,JJf) ! use fine grid topo, note topo diff between coarse/fine grid 
!           topo_sp(i,j)=topoc(IIc,JJc) ! use coarse grid topo, note topo diff between coarse/fine grid 
           latrot_sp(j)=latrot(JJf)
           lonrot_sp(i)=lonrot(IIf)
           topofmc(ista,ijsp)=topo(IIf,JJf)-topoc(IIc,JJc)

           if (dbg) then
             print*, 'source point i,j= ', i,j
             print*, 'source point IIc, JJc= ', IIc, JJc
             print*, 'source point IIf, JJf= ', IIf, JJf
             print*, 'topof,topoc= ',  topo(IIf,JJf), topoc(IIc,JJc)
             print*, 'topof-topoc= ',  topo(IIf,JJf)-topoc(IIc,JJc)
           endif
        ENDDO
        ENDDO
! -----------------------------------------------------------------------
!*** For the set of four source points, call MD to calculate
!*   S and s from each of the four source points
        nobsng=nxy_sp
        allocate(stnsng(nobsng))   
        allocate(xyzsng_rot_8(nobsng,3))   
        allocate(latrotsng(nobsng),lonrotsng(nobsng))   
        allocate(hgtsng(nobsng),STAT=istatus)   
        allocate(qualsng(nobsng),STAT=istatus)
  
        allocate(Sdaughter_3D(nx,ny,nobsng),STAT=istatus)   
        allocate(rg_3D(nx,ny,nobsng),STAT=istatus)        
        allocate(Rnearest_HA(nobsng),STAT=istatus)
        allocate(Znearest_HA(nobsng),STAT=istatus)
        allocate(i_HAx(nobsng),STAT=istatus)
        allocate(j_HAy(nobsng),STAT=istatus)
        allocate(Sdaughter_3Dp(nx,ny,nobsng),STAT=istatus)   
  
!*** treat the four neighbouring coarse grid points as a set of source points
        stnsng='grido'
        iobsng=0
        DO j=1,ny_sp
        DO i=1,nx_sp
           iobsng=iobsng+1
           DO k=1,3
              xyzsng_rot_8(iobsng,k)=xyzsp_rot_8(i,j,k)
           ENDDO
           latrotsng(iobsng)=latrot_sp(j)
           lonrotsng(iobsng)=lonrot_sp(i)
           hgtsng(iobsng)=topo_sp(i,j)
           qualsng(iobsng)=100.          ! -999. for missing incr_sng
        ENDDO
        ENDDO

        runname='test'

        CALL init_para(mothera,motherb,zref1,zref2,avehpbl,runname)
        CALL MD_sph(nobsng,stnsng,lonrotsng,latrotsng,xyzsng_rot_8,hgtsng,  &
                       nx,ny,lonrot,latrot,xyzgr_rot_8,topo,                &
                       Rxy,dx,dy,mothera,motherb,zref1,zref2,avehpbl,       &
                       Sdaughter_3D,rg_3D,i_HAx,j_HAy,                      &
                       Rnearest_HA,Znearest_HA,Sdaughter_3Dp)

        if (dbg) then
          write(6,*) 'ista = ', ista
          write(6,*) 'i_HAx,j_HAy'
          write(6,*) (i_HAx(iobsng),j_HAy(iobsng),iobsng=1,nobsng)
          write(6,*) 'Rnearest_HA, Znearest_HA'
          write(6,*) (Rnearest_HA(iobsng),Znearest_HA(iobsng),iobsng=1,nobsng)
        endif

!!*** to calculate sharing factors and correlation from each of the four 
!!    coarse source points (around the obs station) to the obs station 
        CALL cal_cor(Sdaughter_obstn(ista,:),rg_obstn(ista,:),       &
                     corr_obstn(ista,:),nobsng,Sdaughter_3D,         &
                     rg_3D,nx,ny,latrot,lonrot,xyzgr_rot_8,topo,Rxy, &
                     latstn_rot(ista),lonstn_rot(ista),              &
                     xyzstn_rot_8(ista,:),elevstn(ista),             &
                     mothera,motherb,zref1,zref2,avehpbl)

!*** deallocate intermediate arrays after calculation for each surface
!    station
        deallocate(topo_sp)
        deallocate(xyzsp_rot_8)
        deallocate(latrot_sp)
        deallocate(lonrot_sp)
  
        deallocate(stnsng)
        deallocate(xyzsng_rot_8)
        deallocate(hgtsng)
        deallocate(latrotsng)
        deallocate(lonrotsng)
        deallocate(qualsng)
  
        deallocate(Sdaughter_3D)
        deallocate(rg_3D)
        deallocate(i_HAx)
        deallocate(j_HAy)
        deallocate(Rnearest_HA)
        deallocate(Znearest_HA)
        deallocate(Sdaughter_3Dp)
  
2000   CONTINUE
      DO ista=1,nobs
         print*, ' ista= ', ista
         print*, 'Sdaughter_stn'
         write(6,*) (Sdaughter_obstn(ista,iobsng),iobsng=1,nobsng)
         print*, 'rg_obstn'
         write(6,*) (rg_obstn(ista,iobsng),iobsng=1,nobsng)
         print*, 'corr_obstn'
         write(6,*) (corr_obstn(ista,iobsng),iobsng=1,nobsng)
      ENDDO

      return
      end

