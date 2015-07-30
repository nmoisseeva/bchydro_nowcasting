      program main
!*********************************************************************
! author:	Xingxiu Deng
! date:		January 2007
! Object:	to be used as a main driver 
!*********************************************************************
! input:        indir: directories containing inputs (len=70)
!               topofile_coarse: file name for coarse grid topo file
!               (len=30)
!               topofile_fine: file name for fine grid topo file
!               (len=30)
!               sfcstafile: file name containing the list of all sfc
!                           stations (len=30)
!               maxsta: maximum number of surface stations allowed
!               Rxy: Radius of influence (km)
! output:
!               outdir: directory to store output (len=70)
!               outfstfile: file name of output sharing factor file in
!                           RPN standard format (len=30)
! Note: 1) all the three input files should reside in indir
!       2) the length of characters must be declaired to be the same
!          as indicated above
!*********************************************************************
      implicit none
      character(len=70) indir
      character(len=30) topofile_coarse, topofile_fine
      character(len=30) sfcstafile
      integer maxsta
      character(len=100) topopf_coarse, topopf_fine
      character(len=100) sfcstapf

      character(len=70) outdir
      character(len=30) outfstfile
      character(len=100) outfstpf
      real Rxy

!***  input
      indir='/users/dor/afsg/xin/home/moda/src/inputs/'
      topofile_coarse='lam_mef_coarse.fst'
      topofile_fine='lam_mef_fine.fst'
      sfcstafile='sfcstation.txt'
      maxsta=2000
      Rxy=75.0

      topopf_coarse=trim(indir)//trim(topofile_coarse)
      topopf_fine=trim(indir)//trim(topofile_fine)
      sfcstapf=trim(indir)//trim(sfcstafile)

!***  output
      outdir='/users/dor/afsg/xin/home/moda/src/output/'
      outfstfile='md_sf.fst'
      outfstpf=trim(outdir)//trim(outfstfile)

!*** call the entry subroutine to the sharing factor generator
      call sf_gen(outfstpf,topopf_coarse,topopf_fine,sfcstapf, &
                      maxsta,Rxy,indir)
      STOP
      END
