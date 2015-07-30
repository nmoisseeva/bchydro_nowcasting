       SUBROUTINE getninj_Z(fstZfile,iun,NI,NJ)
!**************************************************************
! author:	Xingxiu Deng
! date:		Oct. 7, 2002 
! Object:	to open fst on Z grid and to get dimension size, 
!               etc. by enquiring positional records ">>", "^^"
!***************************************************************
!      r.compile -src *.f90  -librmn
!--------------------------------------------------------------
!Input variables
!      fstZfile: input FST path and file name
!--------------------------------------------------------------
!Output variables
!      NI: x-dimension size of the input fields
!      NJ: y-dimension size of the input fields
!**************************************************************
      IMPLICIT none
      include 'fst_variables.cdk'

!*---- Defining Mic. Variables
      integer :: i,j,k
      integer :: iun

!*---- Defining input Z grid file
      character (len=100) fstZfile

!*--- output
      INTEGER :: NI,NJ

!-----------------------------------------------------
! debug flag
!-----------------------------------------------------
       logical dbg         ! flag to turn on/off debug
       parameter (dbg = .true.)

!======  Begining of excecutable code

!*--------------------------------------------------
!*     open FST Z grid file to get ni, nj
!*--------------------------------------------------
      print *,"-- Opening FST Z grid File --",trim(fstZfile)
      ier=fnom(iun,trim(fstZfile),'RND',0)

      if (ier.lt.0) then
        print*,'Fatal error while opening the file: ',trim(fstZfile)
        goto 999
      endif

      nrecs = fstouv (iun, 'RND')
      if (nrecs.lt.1) then
         print*,'Error, no record found in file: ',trim(fstZfile)

         goto 999
      endif

!*--------------------------------------------------
!*     Get grid dimension size
!*--------------------------------------------------
      print *,"---- Get grid definition ----"
      keyX = fstinf(iun,nc,nr,nl,-1,' ',-1,-1,-1,' ','>>')
      ni = nc
      ier=fstprm(keyX,dateo,deet,npas,nc,nr,nl,nbits,datyp, &
                   ip1g,ip2g,ip3g,                          &
                   typvar,varname,etiket,                   &
                   grtypg, ig1g,ig2g,ig3g,ig4g,             &
                   swa,lng,dltf,ubc,extra1,extra2,extra3)

      keyY = fstinf(iun,nc,nr,nl,-1,' ',                    &
                   -1,-1,-1,' ','^^')
       nj = nr

      print *,"Output grid ni,nj: ",ni,nj

      IF (dbg) THEN
       print*, 'ip1,ip2,ip3',ip1g,ip2g,ip3g
       print*, 'ig1,ig2,ig3,ig4',ig1g,ig2g,ig3g, ig4g
       print*, 'grtyp ', grtypg
      ENDIF

      if (ni*nj .le. 0) then
         ier = fstfrm(iun)
         ier = fclos(iun)
         print *,"ERROR: grid size = 0"
         print *,"---- ABORT ----"
         goto 999
      endif
 
 999  continue
! Close FST file
      ier = fstfrm(iun)
      ier = fclos(iun)

      RETURN
      END SUBROUTINE getninj_Z



