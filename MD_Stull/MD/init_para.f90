      subroutine init_para(a,b,zref1,zref2,avehpbl,runname)
!***********************************************
!* This SUBROUTINE is to set or read parameters 
!* from namelist.
!*                 Xingxiu Deng, Jan. 17, 2003
!***********************************************
        IMPLICIT none
        REAL :: a,b
        REAL :: zref1, zref2
        REAL :: avehpbl
        CHARACTER (LEN=40) :: runname
        CHARACTER (LEN=80) :: param_file

        namelist /mdpara/a, b, zref1, zref2,avehpbl

!* set default values
        a = 2.
        b = 2.
        zref1 = 750.    ! m, same unit as topo
        zref2 = 750.    ! m 
        avehpbl = 1500.    ! m 
        param_file = 'param_'//trim(runname)
        print*, 'runname is ', trim(runname)

!* read parameters from namelist file
        print*, 'read free parameters from ', trim(param_file),'-end-'
        OPEN(30, file=trim(param_file),action='read', status='old')
        read(30, mdpara, end =100)
        close(30)
        go to 200
100     print*, 'Error reading NAMELIST file - param_input. Default values used'
200     CONTINUE
        print*, 'a = ',a
        print*, 'b = ',b
        print*, 'zref1 = ',zref1
        print*, 'zref2 = ',zref2
        print*, 'avehpbl = ',avehpbl

        return
        end
