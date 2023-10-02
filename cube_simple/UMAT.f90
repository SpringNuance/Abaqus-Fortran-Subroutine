! --------------------
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
     RPL,DDSDDT,DRPLDE,DRPLDT,&
     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
!
     INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME 
      DIMENSION STRESS(NTENS),STATEV(NSTATV), &
     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS), &
     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1), &
     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
!
! -----------------------------------------------------------
!           UMAT FOR ISOTROPIC ELASTICITY
!           CANNOT BE USED FOR PLANE STRESS
! -----------------------------------------------------------
!        ELASTIC PROPERTIES
!
    one = 1.0	
    two = 2.0
    three = 3.0
!
! ===========================================================
! elastic constants 
! ===========================================================
!
    EMOD=PROPS(1) 
    ENU=PROPS(2) 
    EBULK3=EMOD/(one-two*ENU)
    EG2=EMOD/(one+ENU) 
    EG=EG2/two 
    ELAM=(EBULK3-EG2)/three
!
! ===========================================================
! elastic stiffness matrix
! ===========================================================
!
    DO K1=1, NDI
        DO K2=1, NDI
            DDSDDE(K2,K1)=ELAM
        END DO
        DDSDDE(K1,K1) = EG2+ELAM
    END DO
!
    DO K1=NDI+1, NTENS
        DDSDDE(K1,K1)=EG
    END DO
!
! ===========================================================
! Calculate stress
! ===========================================================
!
    DO K1=1, NTENS
        DO K2=1, NTENS
            STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*DSTRAN(K1)
        END DO
    END DO
!
! ===========================================================
!
    RETURN
    END