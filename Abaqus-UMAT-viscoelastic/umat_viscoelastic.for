c ======================================================================
c User Subroutine UMAT for Abaqus viscoelastic material
c By Irfan Habeeb CN (PhD, Technion - IIT)
c using the formulation https://imechanica.org/files/Kelvin-Voigt-Code-Development_0.pdf
c ======================================================================
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
C
      integer k1, k2
      real E, nu, lambda, mu, S(6), D1(6,6), D2(6,6), D3(6,6)

C material properties
      E = props(1)      ! Young's modulus
      nu = props(2)      ! Poisson's ratio
      eta = props(3)    ! viscoelastic coef.

C Lame's parameters
      lambda = E*nu/((1.0d0+nu)*(1.0d0-2.0d0*nu))
      mu = E/(2.0d0*(1.0d0+nu))

c stiffness matrix with viscous effects
      D1 = 0.d0
      D2 = 0.d0
      D3 = 0.d0
      if (dtime .gt. 0.d0) then
        do k1 = 1, ndi
          do k2 = 1, ndi
            D1(k1, k2) = lambda*(1.0d0 + 3.0d0*eta/(E*dtime))
            D2(k1, k2) = lambda
          end do 
          D1(k1, k1) = lambda*(1.0d0 + 3.0d0*eta/(E*dtime)) +
     1      2.0d0*mu*(1.0d0 + eta/(mu*dtime))
          D2(k1, k1) = lambda + 2.0d0*mu
          D3(k1, k1) = -1.0d0
        end do 
c shear stress
        do k1 = ndi+1, ndi+nshr
          D1(k1, k1) = mu*(1.0d0 + eta/(mu*dtime))
          D2(k1, k1) = mu
          D3(k1, k1) = -1.0d0
        end do 
      else 
        do k1 = 1, ndi
          do k2 = 1, ndi
            D1(k1, k2) = lambda
            D2(k1, k2) = lambda
          end do 
          D1(k1, k1) = lambda + 2.0d0*mu
          D2(k1, k1) = lambda + 2.0d0*mu
          D3(k1, k1) = -1.0d0
        end do 
        do k1 = ndi+1, ndi+nshr
          D1(k1, k1) = mu
          D2(k1, k1) = mu
          D3(k1, k1) = -1.0d0
        end do 
      end if 

c Stiffness matrix
      ddsdde(1:ntens, 1:ntens) = D1(1:ntens, 1:ntens)

c stress in the previous step
      S(1:ntens) = stress(1:ntens)

c Stress increment 
      do k1 = 1, ntens
        do k2 = 1, ntens
          stress(k1) = stress(k1) + D1(k1, k2) * dstran(k2) +
     1      D2(k1, k2) * stran(k2) + D3(k1, k2) * S(k2)
        end do 
      end do 
c
      return
      end
