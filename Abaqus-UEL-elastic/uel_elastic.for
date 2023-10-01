c ======================================================================
c User Subroutine UEL for linear elastic material
c All rights of reproduction or distribution in any form are reserved.
c By Irfan Habeeb CN (PhD, Technion - IIT)
c ======================================================================
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,
     3 npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)

      include 'aba_param.inc'

! ******************** state variables in the code *********************
! svars(1) = exx ( at the instant )
! svars(2) = eyy
! svars(3) = ezz (= 0)
! svars(4) = exy
! svars(5) = sxx
! svars(6) = syy
! svars(7) = szz (= 0)
! svars(8) = sxy
! **********************************************************************
      parameter(ngauss=1, nelem=185, nsdv=8)

      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),
     1 svars(nsvars),energy(8),props(*),coords(mcrd,nnode),
     2 u(ndofel),du(mlvarx,*),v(ndofel),a(ndofel),time(2),
     3 params(3),jdltyp(mdload,*),adlmag(mdload,*),
     4 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

       real*8 aintw(3),xii(2,3),xi(2),dNdxi(3,2),dxdxi(2,2),dxidx(2,2),
     * dNdx(3,2),B(3,6),BT(6,3),eps(3),stn(4),stress(3),C(3,3), 
     * ddsdde(4,4),DB(3,6),sig2(3),sdv(nsvars), shape_f(3),str(3), 
     * dstran(4),B_phase(2,3),sig(4)
       real*8 E,nu,we,deter
      
       real*8 ceps(3),auxm3(3,3),B_B(3,3),aux(6),auxm2(6,3),stiff2(3,6), 
     * antn(3,3),F1(3,3),F2(3),F3(3),F4(3), dv(6), dT,
     * epsceps,phase,aka,alc,gc,phase0,epsi,d_dot, lamda, mu
       integer I, J, K, k1, k2, k3

      common/custom/uvars(nelem, 24, ngauss)
      
      E = props(1)
      nu = props(2)
      lamda = E*nu/((1.d0 + nu)*(1.d0 - 2.d0*nu))
      mu = E/(2.d0*(1.d0 + nu))

      do k1 = 1, ndofel
        do k2 = 1, nrhs
          rhs(k1, k2) = 0.d0
        end do 
        do k2 = 1, ndofel
          amatrx(k2, k1) = 0.d0
        end do 
      end do 

      do k1 = 1, 3
        aintw(k1) = 0.5d0
        shape_f(k1) = 0.d0
      end do 

      do k1 = 1, 2
        do k2 = 1, 3
          xii(k1, k2) = 1.d0/3.d0
        end do 
      end do 

      xi(1) = xii(1, 1)
      xi(2) = xii(2, 1)

      dNdxi = 0.d0
      dNdxi(1, 1) = 1.0
      dNdxi(2, 2) = 1.0
      dNdxi(3, 1) = -1.0
      dNdxi(3, 2) = -1.0

      do k1 = 1, 2
        do k2 = 1, 2
          dxdxi(k1, k2) = 0.d0
          do k3 = 1, nnode
            dxdxi(k1, k2) = dxdxi(k1, k2) + coords(k1,k3)*dNdxi(k3, k2)
          end do 
        end do 
      end do 

      deter = 0.d0
      deter = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
      if (deter .lt. 0.d0) then
        print*, 'Negative Jacobian', deter
        call xit
      end if 

      dxidx(1, 1) = dxdxi(2, 2)/deter
      dxidx(2, 2) = dxdxi(1, 1)/deter
      dxidx(1, 2) = -dxdxi(1, 2)/deter
      dxidx(2, 1) = -dxdxi(2, 1)/deter

      do k1 = 1, nnode
        do k2 = 1, 2
          dNdx(k1, k2) = 0.d0
          do k3 = 1, 2
            dNdx(k1, k2) = dNdx(k1, k2) + dNdxi(k1, k3)*dxidx(k3, k2)
          end do 
        end do 
      end do 

      iy = 0
      do k1 = 1, nnode
        ix = iy + 1
        iy = ix + 1
        B(1, ix) = dNdx(k1, 1)
        B(1, iy) = 0.d0
        B(2, ix) = 0.d0
        B(2, iy) = dNdx(k1, 2)
        B(3, ix) = dNdx(k1, 2)
        B(3, iy) = dNdx(k1, 1)
      end do 

      do k1 = 1, 3
        do k2 = 1, 6
          BT(k2, k1) = B(k1, k2)
        end do 
      end do 

      eps = 0.d0
      do k1 = 1, 3
        do k2 = 1, 6
          eps(k1) = eps(k1) + B(k1, k2) * DU(k2, 1)
        end do 
      end do 

      sdv = 0.d0
      do k1 = 1, nsdv
        sdv(k1) = svars(k1)
      end do 

      C = 0.d0
      C(1, 1) = lamda + 2.d0*mu
      C(1, 2) = lamda
      C(2, 1) = C(1, 2)
      C(2, 2) = C(1, 1)
      C(3, 3) = mu

      stress = 0.d0
      do k1 = 1, 3
        do k2 = 1, 3
          stress(k2) = stress(k2) + C(k2, k1)*eps(k1)
        end do 
      end do 

      sig(1) = stress(1)
      sig(2) = stress(2)
      sig(3) = 0.d0
      sig(4) = stress(3)

      stn(1) = eps(1)
      stn(2) = eps(2)
      stn(3) = 0.d0
      stn(4) = eps(3)

      DB = 0.d0
      do k1 = 1, 3
        do k2 = 1, 6
          do k3 = 1, 3
            DB(k1, k2) = DB(k1, k2) + C(k1, k3)*B(k3, k2)
          end do 
        end do 
      end do 

      we = 0.d0
      we = aintw(1) * dabs(deter)

      do k1 = 1, 6
        do k2 = 1, 3
          rhs(k1, 1) = rhs(k1, 1) - we*BT(k1, k2)*stress(k2)
        end do 
      end do 

      do k1 = 1, 6
        do k2 = 1, 6
          do k3 = 1, 3
            amatrx(k1, k2) = amatrx(k1, k2) + we*BT(k1, k3)*DB(k3, k2)
          end do 
        end do 
      end do 

      sdv(1:4) = stn(1:4)             ! strain
      sdv(5:8) = sig(1:4)             ! stress
      svars(1:nsdv) = sdv(1:nsdv)
      uvars(jelem, 1:nsdv, 1) = sdv(1:nsdv)

      return
      end

************************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
C
      include 'aba_param.inc'
C
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      parameter (ngauss = 1, nelem=185, nsdv=8)
      real*8 E, nu, lamda, mu
      common/custom/uvars(nelem, 24, ngauss)

      E = props(1)
      nu = props(2)
      lamda = E*nu/((1.d0 + nu)*(1.d0 - 2.d0*nu))
      mu = E/(2.d0*(1.d0 + nu))

      do k1 = 1, 2
        do k2 = 1, 2
          ddsdde(k2, k1) = lamda
        end do 
        ddsdde(k1, k1) = lamda + 2.d0*mu
      end do
      ddsdde(3,3) = mu

      do k1 = 1, ntens
        do k2 = 1, ntens
          stress(k2) = stress(k2) + ddsdde(k2, k1)*dstran(k1)
        end do 
      end do 

      kelem = noel - nelem
      do k1 = 1, nsdv
        statev(k1) = uvars(kelem, k1, npt)
      end do 

      return
      end 
