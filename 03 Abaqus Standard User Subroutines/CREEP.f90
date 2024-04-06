subroutine CREEP(decra,deswa,statev,serd,ec,esw,p,qtild, &
    temp,dtemp,predef,dpred,time,dtime,cmname,leximp,lend, &
    coords,nstatv,noel,npt,layer,kspt,kstep,kinc)
!
    include 'aba_param.inc'
!
    character*80 cmname
!
    dimension decra(5),deswa(5),statev(*),predef(*),dpred(*), &
              time(3),ec(2),esw(2),coords(*)

    ! user coding to define decra,  deswa

    return
    end

! https://help.3ds.com/2024/English/DSSIMULIA_Established/SIMACAESUBRefMap/simasub-c-creep.htm