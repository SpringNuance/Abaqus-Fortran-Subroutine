      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'ts.inc'

      character(len=1), parameter :: path_sep='/'
      character(len=4), parameter :: mapf_ext='.map'
      character(len=256) :: OUTDIR, JOBNAME, MAPFILE
      character(len=80) :: HEADER

      if (LOP .eq. 0) then
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
              dprop(i_LE11,i) = 0.0
 100      continue
          CALL STDB_ABQERR(-3,'ts.f: nelcohmax to small', 0, 0.0,
     +        '       ')
 101      continue
          nelmap = i-1
          write(*,*) 'Read', nelmap, 'elements'
          CLOSE(UNIT=100)
          RETURN
      endif

      return

 333  continue
      write(*,*) 'ts.f: UEXT error'
      CALL STDB_ABQERR(-3,'ts.f: UEXT error', 0, 0.0, '       ')

      END

      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)

      INCLUDE 'ABA_PARAM.INC'
      include 'ts.inc'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      if (NSTATV .lt. i_LE11) STOP
      if (NFIELD .lt. i_LE11) STOP
      if (CMNAME /= 'COH') STOP

      do 100 i=1, nelmap
          if (kelmap(2,i) == NOEL) GOTO 101
 100  continue
      CALL STDB_ABQERR(-3,'ts.f: unmapped cohesive element %I', NOEL, 0.0,
     +        '       ')
 101  continue

      FIELD(i_LE11) = dprop(i_LE11,i)
      STATEV(i_LE11) = dprop(i_LE11,i)

      RETURN
      END

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
      include 'ts.inc'

      CHARACTER*80 CMNAME,ORNAME
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
      DIMENSION ARRAY(15),JARRAY(15)
      CHARACTER*3 FLGRAY(15)

      if (CMNAME /= 'AISI') STOP

      CALL GETVRM('LE',
     +            ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)

C If error, write comment to .log file
      IF(JRCD.NE.0)THEN
        WRITE(*,*) 'REQUEST ERROR IN UVARM FOR ELEMENT NUMBER ',
     1      NOEL,'INTEGRATION POINT NUMBER ',NPT
        CALL XIT
      ENDIF

c set uvar = LE11
      UVAR(i_LE11) = ARRAY(I_11)

c save in common block

c is NOEL in kelmap(1,*)?
      do 100 i=1, nelmap
          if (kelmap(1,i) == NOEL) GOTO 101
 100  continue

C   no: return
      RETURN

C   yes: save LE11
 101  continue
C   FIXME: check for integration point number
C   FIXME: valid only for CPxxR
      dprop(i_LE11,i) = ARRAY(I_11)
      END
