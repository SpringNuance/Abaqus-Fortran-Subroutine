C*** 2022 June 15
C
C Model Parameters
C      
C kxdim       =   Problem's dimension (2 in 2D, 3 in 3D problems)
C kxnodel     =   Nodes per element (4 in 4-node, 8 in 8-node elements)
C kxelpernode =   Max number of elements belonging to a node (allocate a large enough number) 
C kxelem      =   Largest element label in the model (may skip numbers)
C kxnode      =   Largest node label in the model (may skip numbers)
C kxmatsdv    =   Number of material state variables      
C         
      module ktransfer
C
      implicit none
   
      integer kxdim, kxnodel, kxelpernode, kxelem, kxnode, kxmatsdv
C
      parameter( kxdim       =    2,
     +           kxnodel     =    4,
     +           kxelpernode =   10,
     +           kxelem      = 5068,
     +           kxnode      = 3426,
     +           kxmatsdv    =   11  )
                    
C
      real*8  PELEM(kxelem),PNODAL(kxnode),GRADP(kxelem,kxdim)
      integer IELCONN(kxelem,kxnodel),INODETOEL(kxnode,kxelpernode)
      
C
      save
C
      end module
C
C***********************************************************************
C
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
C
      use ktransfer
      include 'vaba_param.inc'
c
C     Contents of i_Array
      parameter( i_int_nTotalNodes     = 1,
     *           i_int_nTotalElements  = 2,
     *           i_int_kStep           = 3,
     *           i_int_kInc            = 4,
     *           i_int_iStatus         = 5,
     *           i_int_lWriteRestart   = 6  )

C     Possible values for the lOp argument
      parameter( j_int_StartAnalysis    = 0,      
     *           j_int_StartStep        = 1,      
     *           j_int_SetupIncrement   = 2,      
     *           j_int_StartIncrement   = 3,      
     *           j_int_EndIncrement     = 4,      
     *           j_int_EndStep          = 5,      
     *           j_int_EndAnalysis      = 6 )     


C     Possible values for i_Array(i_int_iStatus)
      parameter( j_int_Continue          = 0,      
     *           j_int_TerminateStep     = 1,      
     *           j_int_TerminateAnalysis = 2)      

C     Contents of r_Array
      parameter( i_flt_TotalTime = 1,
     *           i_flt_StepTime  = 2,
     *           i_flt_dTime     = 3 )
C
      dimension i_Array(niArray),r_Array(nrArray)
C        
      ktotalnodes = i_Array(i_int_nTotalNodes)
      ktotalel    = i_Array(i_int_nTotalElements)
      kStepp      = i_Array(i_int_kStep)
      kIncr       = i_Array(i_int_kInc)
      kmaxel      = kxelem 
      kmaxnodes   = kxnode
C
C      
C***  Start of analysis
      if (lOp .eq. 0) then
        PELEM  = 0.D0
        PNODAL = 0.D0
        GRADP  = 0.D0
      end if     
C      
C***  Start of step
      if (lOp .eq. 1) then
C     find which elements a node belongs to
        CALL KNODETOELCON(ktotalel,ktotalnodes,kmaxel,kmaxnodes)
      end if
C
C***  End of increment
      if (lOp .eq. 4) then
C       calculate nodal pressures
        CALL KPNODAL(kmaxel,kmaxnodes)
      end if
C
 1001 FORMAT(1P8E13.5)
 1002 FORMAT(10I8)
      end
C
C***********************************************************************
C
      SUBROUTINE KNODETOELCON(ktotalel,ktotalnodes,kmaxel,kmaxnodes)
C
C*** For each node in the mesh,
C    finds the elements the node belongs to and
C    stores data in array INODETOEL, which is globally available
C
      use ktransfer      
      include 'vaba_param.inc'
C      
      CHARACTER(len=256) OUTDIR, AELREAD, AELWRITE, ANDREAD, ANDWRITE,
     +                   ANELWRITE
C      
C*** Utility subroutine that returns current output directory
C*** (current as opposed to scratch)      
      CALL VGETOUTDIR(OUTDIR,LENOUTDIR)
C
      IELCONN=0
      AELREAD = TRIM(OUTDIR)//'\elements'
      OPEN(UNIT=105,FILE=AELREAD,STATUS='UNKNOWN')
C*** Read element connectivities
      READ(105,*)
      READ(105,*)
      DO i=1,ktotalel
        READ(105,*) IELEM,(IELCONN(IELEM,m),m=1,kxnodel)
      END DO
      CLOSE(105)
C
C*** Find the elements each node belongs to
      INODETOEL = 0
      DO IELEM = 1,kmaxel
        IF (IELCONN(IELEM,1).NE.0) THEN
          DO INODE=1,kxnodel
            KNODE = IELCONN(IELEM,INODE)
            L = INODETOEL(KNODE,1) + 1
            INODETOEL(KNODE,1) = L
            INODETOEL(KNODE,1 + L) = IELEM          
          END DO          
        END IF    
      ENDDO
C      
      END
C
C***********************************************************************
C      
      subroutine vucharlength(
c Read only variables-
     1     nblock, nfieldv, nprops, ncomp, ndim, nnode, nstatev,
     2     kSecPt, kLayer, kIntPt, jElType, jElem,
     3     totalTime, stepTime, dt,
     4     cmname, coordMp, coordNode, direct, T, props,
     5     field, stateOld,
c Write only variables-
     6     charLength )
c
C*** Read current nodal coordinates for each element and
C    store them as state variables.
C    These values are needed for the calculation of the pressure gradient.
C
      use ktransfer      
      include 'vaba_param.inc'
c
      dimension jElType(3), jElem(nblock), coordMp(nblock,ndim), 
     1          coordNode(nblock,nnode,ndim),
     2          direct(nblock,3,3), T(nblock,3,3), props(nprops),
     3          stateOld(nblock,nstatev), charLength(nblock,ncomp), 
     4          field(nblock, nfieldv)
c
      character*80 cmname
C
C*** Read and store current nodal coordinates for the calculation of gradp
C      
      nold = kxmatsdv
      do k = 1, nblock
        do inode=1,kxnodel
          j=(inode-1)*ndim
          do i=1,ndim
            stateOld(k,nold+j+i) = coordNode(k,inode,i)
          end do
        end do
      enddo
c
      return
      end
C
C***********************************************************************
C        
      SUBROUTINE KPNODAL(kmaxel,kmaxnodes)
C
C*** Calculates average nodal pressure values
C
      use ktransfer      
      include 'vaba_param.inc'
C
      DO INODE=1,kmaxnodes
        N=INODETOEL(INODE,1)
        IF (N.NE.0) THEN
          SUM = 0.D0
          DO I=2,1+N
            SUM = SUM + PELEM(INODETOEL(INODE,I))
          END DO
          PNODAL(INODE) = SUM/(1.D0*N)
        END IF
      ENDDO
C
      END
C
C***********************************************************************
C        
      SUBROUTINE KGRADP2D(NOEL,X,Y,PELNODAL)
C
C*** Calculates pressure gradient at (0,0) for 2-D isoparametric
C    4-node quadrilateral element
C
      use ktransfer      
      include 'vaba_param.inc'
      DIMENSION GRADPR(2),X(4),Y(4),PELNODAL(4)
      DIMENSION AJINV(2,2),DN(2,4),BMTRX(2,4)
C
      AJ11 = 0.25D0*( - X(1) + X(2) + X(3) - X(4) )
      AJ12 = 0.25D0*( - Y(1) + Y(2) + Y(3) - Y(4) ) 
      AJ21 = 0.25D0*( - X(1) - X(2) + X(3) + X(4) )
      AJ22 = 0.25D0*( - Y(1) - Y(2) + Y(3) + Y(4) )
      AJ0 = AJ11*AJ22 - AJ21*AJ12
C
      AJINV(1,1) =   AJ22/AJ0
      AJINV(2,2) =   AJ11/AJ0
      AJINV(1,2) = - AJ12/AJ0
      AJINV(2,1) = - AJ21/AJ0
C
      DN(1,1)    = - 0.25D0
      DN(1,2)    =   0.25D0
      DN(1,3)    =   0.25D0
      DN(1,4)    = - 0.25D0
C
      DN(2,1)    = - 0.25D0
      DN(2,2)    = - 0.25D0
      DN(2,3)    =   0.25D0
      DN(2,4)    =   0.25D0
C
      CALL KMULT(AJINV,DN,BMTRX,2,2,4)
      CALL KMULT(BMTRX,PELNODAL,GRADPR,2,4,1)
C
      GRADP(NOEL,1) = GRADPR(1)
      GRADP(NOEL,2) = GRADPR(2)
C
      END
C
C***********************************************************************
C
      SUBROUTINE KGRADP3D(NOEL,X,Y,Z,PELNODAL)
C
C*** Calculates pressure gradient at (0,0,0) for 3-D isoparametric
C    8-node hexahedral element
C
      use ktransfer      
      include 'vaba_param.inc'
      DIMENSION GRADPR(3),X(8),Y(8),Z(8),PELNODAL(8)
      DIMENSION KSI(8),ETA(8),ZETA(8)
      DIMENSION AF(3,3),AFINV(3,3),AFINVT(3,3),DN(3,8),BMTRX(3,8)
C
C*** Nodal coordinates of element in natural space
C
      KSI(1) = -1.D0;  ETA(1) = -1.D0;  ZETA(1) = -1.D0;
      KSI(2) =  1.D0;  ETA(2) = -1.D0;  ZETA(2) = -1.D0;
      KSI(3) =  1.D0;  ETA(3) =  1.D0;  ZETA(3) = -1.D0;
      KSI(4) = -1.D0;  ETA(4) =  1.D0;  ZETA(4) = -1.D0;
C
      KSI(5) = -1.D0;  ETA(5) = -1.D0;  ZETA(5) =  1.D0; 
      KSI(6) =  1.D0;  ETA(6) = -1.D0;  ZETA(6) =  1.D0;
      KSI(7) =  1.D0;  ETA(7) =  1.D0;  ZETA(7) =  1.D0;
      KSI(8) = -1.D0;  ETA(8) =  1.D0;  ZETA(8) =  1.D0;
C
C*** Calculate [F] (Jacobian matrix) at (0,0,0) 
C    
      AF = 0.D0
      DO INODE=1,8
        AF(1,1) = AF(1,1) + KSI(INODE)*X(INODE) 
        AF(2,1) = AF(2,1) + KSI(INODE)*Y(INODE)
        AF(3,1) = AF(3,1) + KSI(INODE)*Z(INODE)
C
        AF(1,2) = AF(1,2) + ETA(INODE)*X(INODE) 
        AF(2,2) = AF(2,2) + ETA(INODE)*Y(INODE)
        AF(3,2) = AF(3,2) + ETA(INODE)*Z(INODE)
C
        AF(1,3) = AF(1,3) + ZETA(INODE)*X(INODE)
        AF(2,3) = AF(2,3) + ZETA(INODE)*Y(INODE)
        AF(3,3) = AF(3,3) + ZETA(INODE)*Z(INODE)
      END DO
      DO I=1,3
      DO J=1,3
        AF(I,J)=AF(I,J)/8.D0
      ENDDO
      ENDDO
C
C*** Calculate [F]^(-T) at (0,0,0)
C
      CALL KINV3X3(AF,AFINV)
      AFINVT = TRANSPOSE(AFINV)
C
C*** Calculate [N'] at (0,0,0)
C
      DO INODE=1,8
        DN(1,INODE) =  KSI(INODE)/8.D0
        DN(2,INODE) =  ETA(INODE)/8.D0
        DN(3,INODE) = ZETA(INODE)/8.D0
      END DO
C
C*** Calculate {gradp}=[F]^(-T).[N'].{PN}  at (0,0,0)
C
      CALL KMULT(AFINVT,DN,BMTRX,3,3,8)
      CALL KMULT(BMTRX,PELNODAL,GRADPR,3,8,1)
C
      GRADP(NOEL,1) = GRADPR(1)
      GRADP(NOEL,2) = GRADPR(2)
      GRADP(NOEL,3) = GRADPR(3)
C
      END
C
C***********************************************************************
C
      SUBROUTINE KMULT(A,B,C,L,M,N)
C
C*** CALCULATES THE MATRIX PRODUCT [C]=[A][B]      
C 
      INCLUDE 'VABA_PARAM.INC'
C
      DIMENSION A(L,M),B(M,N),C(L,N)
C
      DO I=1,L
      DO J=1,N
        AUX = 0.D0
        DO K=1,M
          AUX = AUX + A(I,K)*B(K,J)
        END DO
        C(I,J) = AUX
      END DO
      END DO
C
      END
C
C***********************************************************************
C              
      SUBROUTINE KINV3X3(A,AINV)
C
C*** CALCULATES THE INVERSE OF A 3X3 MATRIX      
C
      INCLUDE 'VABA_PARAM.INC'
C
      DIMENSION A(3,3),AINV(3,3)
C
      DET =  A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3))
     +     - A(1,2)*(A(2,1)*A(3,3) - A(3,1)*A(2,3))
     +     + A(1,3)*(A(2,1)*A(3,2) - A(3,1)*A(2,2))
      ANORM = DSQRT( A(1,1)*A(1,1) + A(1,2)*A(1,2) + A(1,3)*A(1,3)
     +             + A(2,1)*A(2,1) + A(2,2)*A(2,2) + A(2,3)*A(2,3)
     +             + A(3,1)*A(3,1) + A(3,2)*A(3,2) + A(3,3)*A(3,3) )
      TOL = ANORM*1.D-16
C
      IF (DABS(DET)<= TOL) THEN ! Mohsen
        WRITE(*,*) 'TRYING TO INVERT SINGULAR 3X3 MATRIX'
        WRITE(*,*) 'PROGRAM STOPS.'
        CALL XPLB_EXIT
      END IF
C
      AINV(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))/DET
      AINV(1,2) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))/DET
      AINV(1,3) =  (A(1,2)*A(2,3) - A(2,2)*A(1,3))/DET
      AINV(2,1) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))/DET
      AINV(2,2) =  (A(1,1)*A(3,3) - A(3,1)*A(1,3))/DET
      AINV(2,3) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))/DET
      AINV(3,1) =  (A(2,1)*A(3,2) - A(3,1)*A(2,2))/DET
      AINV(3,2) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))/DET
      AINV(3,3) =  (A(1,1)*A(2,2) - A(2,1)*A(1,2))/DET
C
      END
C          
C***********************************************************************
C
      subroutine vumatht (
C Read only (unmodifiable) variables -
     *     nblock, nElem, nIntPt, nLayer, nSectPt, 
     *     ntgrad, nstatev, nfieldv, nprops,  
     *     cmname, stepTime, totalTime, dt,  
     *     coordMp, density, props,  
     *     tempOld, fieldOld, stateOld, enerThermOld, 
     *     tempNew, tempgradNew, fieldNew, 
C Write only (modifiable) variables -
     *     stateNew, fluxNew, enerThermNew, dEnerThDTemp, condEff )
C
      use ktransfer 
      include 'vaba_param.inc'
C
      dimension nElem(nblock)
C
      dimension coordMp(nblock,*), density(nblock), props(nprops),
     *    tempOld(nblock), fieldOld(nblock, nfieldv), 
     *    stateOld(nblock, nstatev), enerThermOld(nblock),
     *    tempNew(nblock), tempgradNew(nblock, ntgrad), 
     *    fieldNew(nblock, nfieldv), stateNew(nblock, nstatev),
     *    fluxNew(nblock, ntgrad), enerThermNew(nblock), 
     *    dEnerThDTemp(nblock,2), condEff(nblock)
C
      character*80 cmname
C
      parameter ( zero = 0.d0, one = 1.d0 )
C
      D     = props(7)
      VH    = props(8)
      AKT   = props(9)
      ALPHA = props(10)
      ANL   = props(11)
      BETA  = props(12)
C
      do k = 1, nblock
C
         noel = nElem(k)
C
         CLT   = tempOld(k)
         CL    = tempNew(k)
         EBAR  = stateOld(k,1)
         RHO   = density(k)
C
C*** Need to account for initial CL and CT          
C
         if (totalTime .eq. zero) then         
           CL0 = CLT
           EBAR0 = zero
           CALL KCT(CT0,THETAL0,THETAT0,EBAR0,CL0,PROPS,NPROPS)
           U0 = (CL0 + CT0)/RHO
           stateOld(k,11) = U0
         end if         
C
C*** Update internal thermal energy and its derivatives
C         
         CALL KCT(CT,THETAL,THETAT,EBAR,CL,PROPS,NPROPS)
         U0 = stateOld(k,11)
         enerThermNew(k) = (CL + CT)/RHO - U0
C        
         dTHLdCL = one/(BETA*ANL)
         AUX = one - THETAL + AKT*THETAL
         dTHTdCL = AKT*dTHLdCL/(AUX*AUX)
         CALL KNT(ANT,DANTDE,EBAR,PROPS,NPROPS)
         dCTdCL = ALPHA*dTHTdCL*ANT
         dEnerThDTemp(k,1) = (one + dCTdCL)/RHO
         dEnerThDTemp(k,2) = zero
C         
C*** Update heat flux vector
C
         aux1 = D*VH
         aux  = aux1*CL
         do i = 1, ntgrad
            fluxNew(k, i) = -D*tempgradNew(k,i) + aux*gradp(noel,i)
         end do
         condEff(k) = D
C
C*** Update the state variables
C       
          stateNew(k,4)  = gradp(noel,1)
          stateNew(k,5)  = gradp(noel,2)
          stateNew(k,6)  = CL
          stateNew(k,7)  = CT
          stateNew(k,8)  = CL + CT
          stateNew(k,11) = U0
C          
      end do
C
      return
      end
C          
C***********************************************************************
C
      subroutine vumat(
c Read only -
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )
c
      include 'vaba_param.inc'
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname
c
      parameter (
     *     i_umt_nblock = 1,
     *     i_umt_npt    = 2,
     *     i_umt_layer  = 3,
     *     i_umt_kspt   = 4,
     *     i_umt_noel   = 5 )
c
c*** Call VUMAT with extra arguments that provide the element numbers
      call  vumatXtrArg(jblock(i_umt_nblock),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(i_umt_noel), jblock(i_umt_npt),
     *     jblock(i_umt_layer), jblock(i_umt_kspt))
c
      end
C
C***********************************************************************
C
      subroutine vumatXtrArg(
c Read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     *     nElement, nMatPoint, nLayer, nSecPoint )
C
      use ktransfer  
      include 'vaba_param.inc'
C
      parameter (zero=0.d0, one=1.d0, two=2.d0 , three=3.d0,
     1 third =1.d0/3.d0, half=0.5d0, op5=1.5d0)
C      
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1     charLength(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)
C
C Documentation of extra arguments:
C  nElement:  Array of internal element numbers
      dimension nElement(nblock)
C  nMatPoint: Integration point number
C  nLayer   : Layer number for composite shells and layered solids
C  nSecPoint: Section point number within the current layer
C
      dimension coornod(kxnodel,kxdim),xc(kxnodel),yc(kxnodel),
     +          zc(kxnodel),pn(kxnodel)
C
      character*80 cmname
C
      iout = 6     ! iout=6 Writes on the .log file
C
C VUMAT with von Mises Plasticity + Hydrogen Effects
C Code for Plane Strain, 3D and Axisymmetric Problems
C
C Material State Variables
C
C state(1)     = ebarp
C state(2)     = yflag
C state(3)     = pbar
C state(4)     = gradpx
C state(5)     = gradpy
C state(6)     = CL
C state(7)     = CT
C state(8)     = CL + CT
C state(9)     = EHKK      
C state(10)    = SMALLC0
C state(11)    = U0
C
C Extra State Variables (needed for the gradp calc.)
C      
C state (12 -- kxmaxsdv + kxnodel*kxdim) = nodal coordinates 
C
C Properties:
C      
C props(1)  = Youngs modulus
C props(2)  = Poisson's ratio
C props(3)  = Initial yield stress
C props(4)  = Hardening exponent
C props(5)  = sig0/e
C props(6)  = Normalizing concentration (C0)
C props(7)  = Lattice diffusion constant (D)
C props(8)  = Partial molar volume of H (VH)
C props(9)  = Equilibrium constant (KT)
C props(10) = H atoms/trap (ALPHA)
C props(11) = Number of host atoms/lattice volume (NL)      
C props(12) = NILS per host atom (BETA)
C props(13) = constant LAMBDA      
C
      e      = props(1)
      xnu    = props(2)
      sig0   = props(3)
      expo   = props(4)      
      e0     = props(5)
      C0     = props(6)
      ANL    = props(11)
      alamda = props(13)
      twomu  = e/(one+xnu)
      thremu = op5*twomu
      alame  = xnu*twomu/(one-two*xnu)
      ak     = e/(three*(one-two*xnu))
C
      if (stepTime.eq.zero) then
        do k=1,nblock
C
          noel = nElement(k)
C            
          trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
          stressNew(k,1) = stressOld(k,1)
     1         +twomu*strainInc(k,1)+ alame*trace
          stressNew(k,2) = stressOld(k,2)
     1         +twomu*strainInc(k,2)+ alame*trace
          stressNew(k,3) = stressOld(k,3)
     1         +twomu*strainInc(k,3)+ alame*trace
          stressNew(k,4) = stressOld(k,4)
     1         +twomu*strainInc(k,4)+ alame*trace
          if(nshr.gt.1) then
            stressNew(k,5)=stressOld(k,5)
     1         +twomu*strainInc(k,5)
            stressNew(k,6)=stressOld(k,6)
     1         +twomu*strainInc(k,6)
          end if
        end do
      else
        do k=1,nblock
C
          noel = nElement(k)
C        
C*** Read nodal coordinates (from vucharlength)
C          
          nold = kxmatsdv
          do inode=1,kxnodel
            j=(inode-1)*kxdim
            do i=1,kxdim
              coornod(inode,i) = stateOld(k,nold+j+i)
            end do
          end do
C
          xc=zero; yc=zero; zc=zero;
C        
          do inode=1,kxnodel
            xc(inode) = coornod(inode,1)
            yc(inode) = coornod(inode,2)
            if (kxdim .gt. 2) then
              zc(inode) = coornod(inode,3)
            end if
          end do
C
C*** Calculate pressure gradient
C          
          do i=1,kxnodel
            pn(i) = PNODAL(IELCONN(noel,i))
          end do
          if (kxdim .eq. 2) then
            call KGRADP2D(noel,xc,yc,pn)
          else 
            call KGRADP3D(noel,xc,yc,zc,pn)
          end if
C          
C*** Form SMALLC0
C          
          if (totalTime .eq. zero) then
            CLT = tempOld(k)
            CALL KCT(CTT,THETAL,THETAT,EBART,CLT,PROPS,NPROPS)
            SMALLC0 = (CLT + CTT)/ANL
            stateOld(k,7)  = CTT
            stateOld(k,10) = SMALLC0
          end if                  
C 
          peeqOld = stateOld(k,1)
          CTT     = stateOld(k,7)
          EHKKT   = stateOld(k,9)
          SMALLC0 = stateOld(k,10)
          CLT     = tempOld(k)
          CL      = tempNew(k)  
C        
C*** Elastic Predictor Calculation
C            
          call kycurve(yieldOld,hard,peeqOld,CL,expo,e0,sig0,props,
     +                 nprops)
          trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
          s11 = stressOld(k,1) + twomu*strainInc(k,1)  + alame*trace    
          s22 = stressOld(k,2) + twomu*strainInc(k,2)  + alame*trace
          s33 = stressOld(k,3) + twomu* strainInc(k,3) + alame*trace
          s12 = stressOld(k,4) + twomu* strainInc(k,4)
          if (nshr.gt.1) then
            s13 = stressOld(k,5) + twomu* strainInc(k,5)
            s23 = stressOld(k,6) + twomu* strainInc(k,6)
          end if
C
          smean = third*(s11+s22+s33)
          s11 = s11-smean
          s22 = s22-smean
          s33 = s33-smean
          if (nshr.eq.1) then
            vmises = sqrt (op5*(s11*s11 + s22*s22 + s33*
     1                    s33 + two*s12*s12))
          else
            vmises = sqrt(op5*(s11*s11 + s22*s22 + s33*
     1                    s33 + two*s12*s12 +two*s13*
     2                    s13 + two*s23*s23))
          end if
C
          sigdif = vmises-yieldOld
          facyld = zero
          yflag  = zero
          phydro = zero
          if (sigdif.gt.zero) then
            facyld = one
            yflag  = one
          end if
          deqps = facyld*sigdif/(thremu + hard)
C
C*** Update stresses
C
          yieldNew = yieldOld + hard*deqps
          peeqNew  = peeqOld + deqps
          factor   = yieldNew/(yieldNew + thremu*deqps)
          if (sigdif.gt.zero) then
            CALL KSMALLC(SMALLC,dcdEP,dcdCL,peeqNew,CL,props,nprops)
            aux = one + alamda*(SMALLC-SMALLC0)/three
            EHKK = three*DLOG(AUX)
            phydro = ak*(EHKK - EHKKT)
          end if
          pNew = smean - phydro
          stressNew(k,1) = s11*factor + pNew
          stressNew(k,2) = s22*factor + pNew
          stressNew(k,3) = s33*factor + pNew
          stressNew(k,4) = s12*factor
          if (nshr.gt.1) then
            stressNew(k,5) = s13*factor
            stressNew(k,6) = s23*factor
          end if
C     
C*** Calculate smooth pressure field pbar
C
          smeanNew=third*(stressNew(k,1)+stressNew(k,2)+stressNew(k,3))
          PELEM(noel) = smeanNew
C      
          sum    = zero
          annode = zero
          do i=1,kxnodel
            pn(i) = PNODAL(IELCONN(noel,i))
            sum   = sum + pn(i)
            annode = annode + one
          end do
          pbar = sum/annode        
C
C*** Update the state variables
C
          stateNew(k,1) = peeqNew
          stateNew(k,2) = yflag
          stateNew(k,3) = -pbar
          stateNew(k,9) = EHKK
C
C*** Update the specific internal energy
C
          if (nshr.eq.1 ) then
            stressPower = half*(
     1      (stressOld(k,1) + stressNew(k,1))*strainInc(k,1)+
     2      (stressOld(k,2) + stressNew(k,2))*strainInc(k,2)+
     3      (stressOld(k,3) + stressNew(k,3))*strainInc(k,3))+
     4      (stressOld(k,4) + stressNew(k,4))*strainInc(k,4)
          else
            stressPower=half*(
     1      (stressOld(k,1) + stressNew(k,1))*strainInc(k,1)+
     2      (stressOld(k,2) + stressNew(k,2))*strainInc(k,2)+
     3      (stressOld(k,3) + stressNew(k,3))*strainInc(k,3))+
     4      (stressOld(k,4) + stressNew(k,4))*strainInc(k,4)+
     5      (stressOld(k,5) + stressNew(k,5))*strainInc(k,5)+
     6      (stressOld(k,6) + stressNew(k,6))*strainInc(k,6) 
          end if
          enerInternNew(k)= enerInternOld(k) + stressPower/density(k)
C
C*** Update the dissipated inelastic specific energy 
C
          plasticWorkInc  = half*(yieldOld+yieldNew)*deqps
          enerInelasNew(k)= enerInelasOld(k)+plasticWorkInc/density(k)
C
        end do
        end if
C
      return
      end
C
C***********************************************************************
C
      SUBROUTINE kycurve(yield,hard,ebar,cl,expo,e0,sig0,props,nprops)
C
C*** Calculates yield stress and hardening modulus of the material
C*** Isotropic power-law hardening with n=expo
C 
      include 'vaba_param.inc'
      dimension props(nprops)
C
C Power law hardening
C
C
      CALL KSMALLC(SMALLC,dcdEP,dcdCL,EBAR,CL,PROPS,NPROPS)
C
      XI   = 0.1D0
      f    = 1.D0 - (1.D0-XI)*SMALLC
      dfdc = -(1.D0-XI)
      f    = 1.D0
      dfdc = 0.D0
C            
      if (expo.gt.50.d0) then
        yield = f*sig0
        hard  = dfdc*dcdEP*sig0
      else
        aux   = sig0*(1.d0+ebar/e0)**(1.d0/expo)
        yield = f*aux
        hard  = dfdc*dcdEP*aux + yield/(expo*(e0+ebar))
      end if
C      
      END
C
C***********************************************************************      
C
      SUBROUTINE KCT(CT,THETAL,THETAT,EBAR,CL,PROPS,NPROPS)
C
C*** computes C_T, theta_L, theta_T
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS)
C
C     read constants
      AKT   = PROPS(9)
      ALPHA = PROPS(10)
      ANL   = PROPS(11)
      BETA  = PROPS(12)
C
      THETAL = CL/(BETA*ANL)
      THETAT = AKT*THETAL/(1.D0 + AKT*THETAL - THETAL)
      CALL KNT(ANT,DANTDE,EBAR,PROPS,NPROPS)
      CT     = ALPHA*THETAT*ANT
C
      END
C
C***********************************************************************      
C
      SUBROUTINE KNT(ANT,DANTDE,EBAR,PROPS,NPROPS)
C
C*** computes N_T and its derivative wrt ebar^p
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS)
C      
C*** TAHA & SOFRONIS FITTING (2001)
C
      C0 = PROPS(6)
      AA = 23.26D0; BB = 2.33D0; CC = 5.5D0
      AUX    = DEXP(-CC*EBAR)
      ANT    = ( 10**(AA - BB*AUX) )/C0
      DANTDE = BB*CC*DLOG(10.D0)*AUX*ANT
C
      END
C
C***********************************************************************      
C
      SUBROUTINE KSMALLC(SMALLC,dcdEP,dcdCL,EBAR,CL,PROPS,NPROPS)
C
C*** computes c = (C_L + C_T)/N_L and
C    its derivatives wrt ebar^p and C_L
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS)
C
      TINY=1.D-14
C      
C     read constants
      AKT   = PROPS(9)
      ALPHA = PROPS(10)
      ANL   = PROPS(11)
      BETA  = PROPS(12)
C      
      CALL KCT(CT,THETAL,THETAT,EBAR,CL,PROPS,NPROPS)
      SMALLC = (CL + CT)/ANL
C
      IF (SMALLC.LE.TINY) THEN
        dcdEP  = 0.D0
        dcdCL  = 0.D0
        GOTO 1010
      END IF      
C
C*** find dcdEP
      CALL KNT(ANT,DANTDE,EBAR,PROPS,NPROPS)
      dcdEP = ALPHA*THETAT*DANTDE/ANL    
C
C*** find dcdCL
      dTHLdCL = 1.D0/(BETA*ANL)
      dTHTdCL = THETAT*dTHLdCL/(THETAL*(1.D0 + AKT*THETAL - THETAL))
      dCTdCL  = ALPHA*dTHTdCL*ANT
      dcdCL   = (1.D0 + dCTdCL)/ANL
C
 1010 CONTINUE      
      END
C
C***********************************************************************
C
      subroutine vdisp(
c Read only variables -
     *   nblock, nDof, nCoord, kstep, kinc, 
     *   stepTime, totalTime, dtNext, dt, 
     *   cbname, jBCType, jDof, jNodeUid, amp,
     *   coordNp, u, v, a, rf, rmass, rotaryI, 
c Write only variable -
     *   rval )
c
      include 'vaba_param.inc'
      parameter( zero = 0.d0, half = 0.5d0, one = 1.d0 )
c
      character*80 cbname
      dimension jDof(nDof), jNodeUid(nblock), 
     *   amp(nblock), coordNp(nCoord,nblock), 
     *   u(nDof,nblock), v(nDof,nblock), a(nDof,nblock), 
     *   rf(nDof,nblock), rmass(nblock), 
     *   rotaryI(3,3,nblock), rval(nDof,nblock)
C      
C*** Impose displacements
C
      if ( jBCType .NE. 0. OR. stepTime .LT. zero ) RETURN
C
      E    = 828.D0
      ANU  = 0.3D0
      G    = E/(2.D0*(1.D0 + ANU))
      AUXN = 3.D0 - 4.D0*ANU
      PI   = 4.D0*DATAN(1.D0)
      TOL  = -PI/180.D0
C      
C*** SIF = MODE I STRESS INTENSITY FACTOR K_I
C*** Set SIFMAX equal to final normalized load
C*** Set TMAX equal to the total analysis time
C      
      SIFMAX = 40.115d0;  TMAX = 32.5d0
      SIF = SIFMAX*totalTime/TMAX
      AUX1 = 0.5D0*SIF/(G*DSQRT(2.D0*PI))
c
      do k = 1, nblock
        X1  = coordNp(1,k) - u(1,k)
        X2  = coordNp(2,k) - u(2,k)
        R   = DSQRT(X1*X1 + X2*X2)
        TH  = DATAN2(X2,X1)
        IF (TH.LT.TOL) TH = TH + PI
        AUX = AUX1*DSQRT(R)*(AUXN - DCOS(TH))
        rval(1,k) = AUX*DCOS(0.5D0*TH)
        rval(2,k) = AUX*DSIN(0.5D0*TH)
      ENDDO      
c
      return       
      end
C
C***********************************************************************
C
 


