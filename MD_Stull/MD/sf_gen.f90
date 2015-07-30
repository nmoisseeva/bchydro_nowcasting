      SUBROUTINE sf_gen(outfstpf,topopf_coarse,topopf_fine,sfcstapf, &
                            maxsta,Rxy,indir)
!*********************************************************************
! author:	Xingxiu Deng
! date:		January 2007
! Object:	entry to sharing factor generator
!*********************************************************************
! input:        topopf_coarse: path/file name for coarse grid topo file
!               topopf_fine: path/file name for fine grid topo file
!               sfcstapf: path/file name containing the list of all sfc
!                           stations
!               Rxy: Radius of influence (km)
!               indir: directory where holds inputs
! output:       outfstpf: path/file name of output sharing factor file in
!                         RPN standard format
!*********************************************************************
      implicit none
!*** output
      character(len=100) outfstpf

!*** input
      character(len=100) topopf_coarse, topopf_fine
      character(len=100) sfcstapf
      integer maxsta
      real Rxy
      character(len=70) indir

      print*, '---- Welcome to sharing factor generator ----------'

      CALL get_inputs(outfstpf,topopf_coarse,topopf_fine,sfcstapf, &
                      maxsta,Rxy,indir)

      RETURN
      end SUBROUTINE sf_gen

