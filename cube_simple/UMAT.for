SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C
C -----------------------------------------------------------
C           UMAT FOR ISOTROPIC ELASTICITY
C           CANNOT BE USED FOR PLANE STRESS
C -----------------------------------------------------------
C        ELASTIC PROPERTIES
C
    one = 1.0	
    two = 2.0
    three = 3.0
C
C ===========================================================
C elastic constants 
C ===========================================================
C
    EMOD = PROPS(1) 
    ENU = PROPS(2) 
    EBULK3 = EMOD/(one-two*ENU)
    EG2 = EMOD/(one+ENU) 
    EG = EG2/two 
    ELAM (EBULK3 - EG2)/three
C
C ===========================================================
C elastic stiffness matrix
C ===========================================================
C
    DO I = 1,NDI
        DO K2 = 1,NDI
            DDSDDE(K1,K2) = ELAM
        END DO
        DDSDDE(K1,K1) = EG2 + ELAM
    END DO
C
    DO K1 = NDI+1,NTENS
        DDSDDE(K1,K1) = EG
    END DO
C
C ===========================================================
C Calculate stress
C ===========================================================
C
    DO K1 = 1,NTENS
        DO K2 = 1,NTENS
            STRESS(K2) = STRESS(K2) + DDSDDE(K2,K1)*DSTRAN(K1)
        END DO
    END DO
C
C ===========================================================
C
    RETURN
    END