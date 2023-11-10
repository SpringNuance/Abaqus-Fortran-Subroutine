      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

      include 'ABA_PARAM.INC'
      include 'hydra.inc'

      DIMENSION TIME(2)

      character(len=1), parameter :: path_sep='/'
      character(len=4), parameter :: jacf_ext='.jac'
      character(len=256) :: OUTDIR, JOBNAME, JACFILE
      character(len=80) :: HEADER

      write(*,1001) 'UEXTDB', lop, TIME, DTIME, KSTEP, KINC
 1001 FORMAT(' ',A,I3,3F7.3,2I4)

      if (LOP .eq. 0) then
          call GETJOBNAME(JOBNAME, LENJOBNAME)
          call GETOUTDIR(OUTDIR, LENOUTDIR)
          write(*,*) 'UEXT ','Start of ',JOBNAME(1:LENJOBNAME)
          write(*,*) 'UEXT ','Start of analysis in ',OUTDIR(1:LENOUTDIR)

          JACFILE=
     +    OUTDIR(1:LENOUTDIR)//path_sep//JOBNAME(1:LENJOBNAME)//jacf_ext

          open(unit=100,status='OLD',err=333, file=JACFILE)
          read(100,'(A80)') HEADER
          do i=1,nelhmax
              read(100,'(A80)',END=101) HEADER
              read(100,*,ERR=333,END=333) kelhlab(i)
              read(100,*,ERR=333,END=333) (kelhinc(j,i), j=1,NELHINC)
              read(100,*,ERR=333,END=333)
     +                ((delhjac(k,j,1,i),j=1,NELHINC),k=1,ND)
          end do
          call STDB_ABQERR(-3,'hydra: nelhmax too small',
     +                     0, 0.0, '        ')
 101      continue
          nelh = i-1
          write(*,*) 'UEXTDB Read', nelh, 'elh elements'
          close(unit=100)

C init press gradient
          do i=1,ND
              do j=1,NELHGPT
                  do k=1,nelh
                      delhgradph(i,j,k) = 0.0
                  end do
              end do
          end do

          return
      endif

      return

 333  continue
      write(*,*) 'hydra: 333 error'
      call STDB_ABQERR(-3,'hydra: 333 error',
     +                 0, 0.0, '        ')

      end

      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)

      INCLUDE 'ABA_PARAM.INC'
      include 'hydra.inc'

      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))

      dimension kndhlab(nndhmax)
      dimension dpressh(nndhmax)
      dimension kincel(NELHINC)

      write(*,*) 'URDFIL'

C
C find current increment in .fil file
C
      CALL POSFIL(KSTEP,KINC,ARRAY,JRCD)
      if (ARRAY(3) .ne. TIME(1)) goto 334

      i = 0
      do while (.true.)
         CALL DBFILE(0,ARRAY,JRCD)
         IF (JRCD .NE. 0) GO TO 110
         KEY=JRRAY(1,2)
         if (KEY .eq. 1) then
            i = i + 1
            kndhlab(i) = jrray(1,3)
            if (jrray(1,6).ne.4) goto 334
         else if (KEY.EQ.12) then
            dpressh(i) = -array(5)
         end if
      end do
 110  CONTINUE
      nnd = i

      do i=1,nelh
          do j=1,NELHINC
              l = kelhinc(j,i)
              do k=1,nnd
                  if (kndhlab(k) .eq. l) then
                     kincel(j) = k
                     goto 220
                  endif
              end do
 220          continue
          end do
          do j=1,NELHGPT
              do k=1,ND
                  dprod = 0.0
                  do l=1,NELHINC
                      dprod=dprod+delhjac(k,l,j,i)*dpressh(kincel(l))
                  end do
                  delhgradph(k,j,i) = dprod
              end do
          end do
      end do
      RETURN

 334  continue
      write(*,*) 'hydra: 334 error'
      call STDB_ABQERR(-3,'hydra: 334 error',
     +                 0, 0.0, '        ')

      END

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
      include 'hydra.inc'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      do i=1,nelh
          if (kelhlab(i) .eq. NOEL) then
              do j=1,ND
                  UVAR(j) = delhgradph(j,NPT,i)
              end do
          endif
      end do

      END
