      
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'conc.inc'

      character(len=1), parameter :: path_sep='/'
      character(len=4), parameter :: mapf_ext='.map'
      character(len=256) :: OUTDIR, JOBNAME, MAPFILE
      character(len=80) :: HEADER

      if ( LOP .eq. 0) then
          CALL GETJOBNAME(JOBNAME, LENJOBNAME)
          CALL GETOUTDIR(OUTDIR, LENOUTDIR)
          write(*,*) 'UEXT ','Start of ',JOBNAME(1:LENJOBNAME)
          write(*,*) 'UEXT ','Start of analysis in ',OUTDIR(1:LENOUTDIR)

          MAPFILE=
     +    OUTDIR(1:LENOUTDIR)//path_sep//JOBNAME(1:LENJOBNAME)//mapf_ext

          OPEN(UNIT=100,STATUS='OLD',ERR=333,FILE=MAPFILE)
          READ(100,'(A80)') HEADER
          do 100 i=1,nelcohmax
              read(100,*,ERR=333,END=101) kelmap(1,i), kelmap(2,i)
              dprop(N_DEC,i) = 1.0
 100      continue
          CALL STDB_ABQERR(-3,'conc.f: nelcohmax to small', 0, 0.0, 
     +        '       ')
 101      continue
          nelmap = i-1
          write(*,*) 'Read', nelmap, 'elements'
          CLOSE(UNIT=100)
          RETURN
      endif

      return

 333  continue
      write(*,*) 'conc.f: UEXT error'
      CALL STDB_ABQERR(-3,'conc.f: UEXT error', 0, 0.0, '       ')

      END

      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)

      INCLUDE 'ABA_PARAM.INC'
      include 'conc.inc'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     if (CMNAME /= 'COH') RETURN

      do 100 i=1, nelmap
          if (kelmap(2,i) == NOEL) GOTO 101
 100  continue
      CALL STDB_ABQERR(-3,'conc.f: unmapped cohesive element %I', NOEL,
     + 0.0, '       ')
 101  continue
       
      FIELD(N_DEC) = dprop(N_DEC,i)
      STATEV(N_DEC) = dprop(N_DEC,i)
C     write(*,*) 'USDFLD', NOEL, NPT, NFIELD, FIELD(N_DEC)

      RETURN
      END

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
      include 'conc.inc'

      CHARACTER*80 CMNAME,ORNAME
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
      DIMENSION ARRAY(15),JARRAY(15)
      CHARACTER*3 FLGRAY(15)

      JERROR = 0

      CALL GETVRM('PE',ARRAY,JARRAY,FLGRAY,JRCD,
     + JMAC,JMATYP,MATLAYO,LACCFLA)
      JERROR = JERROR + JRCD
      PEEQ = ARRAY(7)

      CALL GETVRM('FV',ARRAY,JARRAY,FLGRAY,JRCD,
     + JMAC,JMATYP,MATLAYO,LACCFLA)
      JERROR = JERROR + JRCD
      F_NNC = ARRAY(N_NNC)

C material constant
      Cl = F_NNC*0.071
C (9)
      Ct = (49.0*PEEQ + 0.1)*Cl
      Ctot = Cl + Ct
C (14)
      theta = Ctot / (Ctot + 4.51)
      dec = 1 - 1.0467*theta + 0.1687*theta*theta

      UVAR(1) = Ct
      UVAR(2) = Ctot
      UVAR(3) = theta
      UVAR(4) = dec

C If error, write comment to .log file
      IF(JERROR.NE.0)THEN
        WRITE(*,*) 'REQUEST ERROR IN UVARM FOR ELEMENT NUMBER ',
     1      NOEL,'INTEGRATION POINT NUMBER ',NPT
        CALL XIT
      ENDIF

C save in common block

C is NOEL in kelmap(1,*)?
      do 100 i=1, nelmap
          if (kelmap(1,i) == NOEL) GOTO 101
 100  continue

C   no: return
      RETURN

C   yes: save dec
 101  continue
      dprop(N_DEC,i) = dec
      END
