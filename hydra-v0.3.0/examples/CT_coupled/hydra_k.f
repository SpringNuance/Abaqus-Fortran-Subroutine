      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

      include 'ABA_PARAM.INC'
      include 'hydra.inc'

      DIMENSION TIME(2)
      character(len=1), parameter :: path_sep='/'
      character(len=4), parameter :: jacf_ext='.jac'
      character(len=4), parameter :: mapf_ext='.map'
      character(len=256) :: OUTDIR, JOBNAME, MAPFILE, JACFILE
      character(len=80) :: HEADER

      if (LOP .eq. 0) then
          call GETJOBNAME(JOBNAME, LENJOBNAME)
          call GETOUTDIR(OUTDIR, LENOUTDIR)
          write(*,*) 'UEXT ','Start of ',JOBNAME(1:LENJOBNAME)
          write(*,*) 'UEXT ','Start of analysis in ',OUTDIR(1:LENOUTDIR)

C read jacobian
          JACFILE =
     +    OUTDIR(1:LENOUTDIR)//path_sep//JOBNAME(1:LENJOBNAME)//jacf_ext
          open(unit=100,status='OLD',err=333,file=JACFILE)
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
          close(unit=100)
          nelh = i-1
          write(*,*) ' UEXTDB elh elements:', nelh

C read map.txt
          MAPFILE=
     +    OUTDIR(1:LENOUTDIR)//path_sep//JOBNAME(1:LENJOBNAME)//mapf_ext

          open(unit=100,status='OLD',err=333,file=MAPFILE)

          read(100,'(A80)') HEADER
          do i=1,nelcohmax
              read(100,*,ERR=333,END=201) kelmap(1,i), kelmap(2,i)
          end do
          CALL STDB_ABQERR(-3,'hydra: nelcohmax too small', 0, 0.0,
     +        '       ')
 201      continue
          nelcoh = i-1
          write(*,*) ' UEXTDB coh elements:', nelcoh
          close(unit=100)

C init press gradient
          do i=1,ND
              do j=1,NELHGPT
                  do k=1,nelh
                      delhgradph(i,j,k) = 0.0
                  end do
              end do
          end do
C init damage
          do k=1,nelcoh
              dcohdam(k) = 1.0
          end do

      else if (LOP .eq. 3) then
C         end of analysis
      else if (LOP .eq. 4) then
C         restart analysis, unsopported
          STOP 'restart unsupported'
      endif

      return

 333  continue
      call STDB_ABQERR(-3,'hydra: 333 error', 0, 0.0, '        ')

      end

      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)

      INCLUDE 'ABA_PARAM.INC'
      include 'hydra.inc'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      if (CMNAME .EQ. 'COH') then

          do i=1, nelcoh
              if (kelmap(2,i) == NOEL) GOTO 101
          end do
          CALL STDB_ABQERR(-3,'hydra: unmapped cohesive element %I',
     +         NOEL, 0.0,
     +        '       ')
 101      continue

          if (NFIELD .lt. N_DAM) then
              CALL STDB_ABQERR(-3,'hydra: USDFLD NFIELD: %I',
     +                         NOEL, 0.0, '       ')
          endif
          if (NSTATV .ne. N_coh_SDV) then
              CALL STDB_ABQERR(-2,'hydra: USDFLD NSTATV: %I',
     +                         NSTATV, 0.0, '       ')
              CALL STDB_ABQERR(-3,'hydra: USDFLD NOEL: %I',
     +                         NOEL, 0.0, '       ')
          endif
          FIELD(N_DAM) = dcohdam(i)
          STATEV(N_DAM) = dcohdam(i)

      else if (CMNAME .EQ. 'AISI') then
C          error checks
           if (NFIELD .NE. NFIELD_AISI) then
              CALL STDB_ABQERR(-3,'hydra: USDFLD AISI NFIELD: %I',
     +                         NFIELD, 0.0, '       ')
           else if (NSTATV .NE. NSTATV_AISI) then
              CALL STDB_ABQERR(-3,'hydra: USDFLD AISI NSTATV: %I',
     +                         NSTATV, 0.0, '       ')
           endif
           FIELD(IPHI_AISI) = STATEV(I_phi)
      else
C          unknown material
           CALL STDB_ABQERR(-3,'hydra: USDFLD CMNAME: %S',
     +                      NOEL, 0.0, CMNAME)
      endif

      RETURN
      END

c urdfil documented in abq usr sub man, 1.1.48
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)

      INCLUDE 'ABA_PARAM.INC'
      include 'hydra.inc'

      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))

      dimension kndhlab(nndhmax)
      dimension dpressh(nndhmax)
      dimension kincel(NELHINC)

      write(*,*) '+URDFIL ', KSTEP, KINC
      write(*,*) ' URDFIL ', TIME

C
C find current increment in .fil file
C

C posfil documented in abq user manual 5.1.4
C record structure see abq user manual 5.1.2
      CALL POSFIL(KSTEP,KINC,ARRAY,JRCD)

c     check end of file condition on .fil
      if (JRCD .ne. 0) goto 334
c     check increment header record (2000)
      if (JRRAY(1,2) .ne. 2000) goto 335
      if (JRRAY(1,8) .ne. KSTEP) goto 335
      if (JRRAY(1,9) .ne. KINC) goto 335

c KEY <- RECORD KEY
c array(i+2)   <- real attribute_i
c jrray(1,i+2) <- integer attribute_i
      i = 0
      do while (.true.)
         CALL DBFILE(0,ARRAY,JRCD)
c        check if end of file
         IF (JRCD .NE. 0) GO TO 110
         KEY=JRRAY(1,2)
         if (KEY .eq. 1) then
            i = i + 1
c save node label
            kndhlab(i) = jrray(1,3)
c attribute_4 == 4: check that results are nodal averaged values
            if (jrray(1,6).ne.4) goto 336
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

c  overwrite FIL data in next increment (KINC+1) of current step (KSTEP)
c  prevents FIL file from growing and problems in finding current increment
c  for very large files
      LOVRWRT = 1
      RETURN

 334  continue
      call STDB_ABQERR(-3,'hydra: 334 error', 0, 0.0, '        ')
      return

 335  continue
      call STDB_ABQERR(-3,'hydra: 335 error', 0, 0.0, '        ')
      return

 336  continue
      call STDB_ABQERR(-3,'hydra: 336 error', 0, 0.0, '        ')
      return

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
      
      CALL GETVRM('PE',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)

      do i=1,nelh
          if (kelhlab(i) .eq. NOEL) then
              goto 101
          endif
      enddo
      call STDB_ABQERR(-3,'hydra: NOEL error', 0, 0.0, '        ')
 101  continue

C save hydro stress gradient
      do j=1,ND
          UVAR(j) = delhgradph(j,NPT,i)
      end do

      if (KINC .eq. 0) then
C save PEEQ at increment start
c fixme: what to do with delhpeeq0 at KINC == 0 ????
          delhpeeq0(NPT,i) = ARRAY(7)
          delhpeeqC(NPT,i) = ARRAY(7)
          delhpeeqCtime(NPT,i) = TIME(2)
      else
          if ( TIME(2) .gt. delhpeeqCtime(NPT,i) ) then
C save current PEEQ
              delhpeeq0(NPT,i) = delhpeeqC(NPT,i)
              delhpeeqC(NPT,i) = ARRAY(7)
              delhpeeqCtime(NPT,i) = TIME(2)
          endif
      endif

      END

      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
      include 'hydra.inc'

      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(1),COORDS(3)

      parameter (RT=2462.)
      parameter (aKt=EXP(60000./RT))
      parameter (aNl=5.1D20, Vh=2.D3)

      if ( NTGRD .ne. ND ) then
          call STDB_ABQERR(-3,'hydra: ND error', 0, 0.0, '        ')
      endif

      do i=1,nelh
          if (kelhlab(i) .eq. NOEL) then
              goto 101
          endif
      enddo
      call STDB_ABQERR(-3,'hydra: NOEL error', 0, 0.0, '        ')
 101  continue

      D = PROPS(1)
      Cl = TEMP
      PEEQ = delhpeeqC(NPT,i)
      DPEEQ = (delhpeeqC(NPT,i)-delhpeeq0(NPT,i))

      aNt = 1.D-9*(10**(23.26-2.33*EXP(-5.5*PEEQ)))
      dNtdep = 29.5*EXP(-5.5*PEEQ)*aNt
      DUDT = 1 + (aNt*aKt*aNl)/(aKt*Cl+aNl)**2
      DU = DUDT * DTEMP + ((aKt*Cl)/(aKt*Cl+aNl))*dNtdep*DPEEQ
      U = U+DU
      do k=1, NTGRD
         FLUX(k) = -D*DTEMDX(k)+(D*Cl*((Vh*1.D-3)/RT)*delhgradph(k,NPT,i))
         DFDT(k) = D*((Vh*1.D-3)/RT)*(delhgradph(k,NPT,i))
         DFDG(k,k) = -D
      end do

      Ct = (aNt*aKt*Cl)/(aKt*Cl+aNl)
      C = (Cl + Ct) * 2.11D-16
c Concentration is in ppm
      theta = C / (C + 5.11)
      dec = 1 - 1.0467*theta + 0.1687*theta*theta

C
C save state variables
C

      if (NSTATV .ne. NSTATV_AISI) then
          CALL STDB_ABQERR(-2,'hydra: UMATHT NSTATV: %I',
     +                     NSTATV, 0.0, '       ')
          CALL STDB_ABQERR(-3,'hydra: UMATHT NOEL: %I',
     +                     NOEL, 0.0, '       ')
      endif
C debug
      STATEV(I_aNt) = aNt
      STATEV(I_Cl) = Cl
      STATEV(I_PEEQ) = PEEQ
      STATEV(I_Ct) = Ct
      STATEV(I_C) = C
      STATEV(I_theta) = theta
      STATEV(I_dec) = dec
C hardening parameter
      STATEV(I_phi) = (Cl/aNl + Ct/aNt)/2.

C is NOEL in kelmap(1,*)?
      do i=1, nelcoh
          if (kelmap(1,i) == NOEL) GOTO 201
      end do

C   no: return
      return

C   yes: save dec and return
 201  continue
      dcohdam(i) = dec
      return

      END
