c ======================================================================
c User Subroutine UMAT for Abaqus linear elastic material
c By Irfan Habeeb CN (PhD, Technion - IIT)
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
      integer i, j
      real Y, n, lambda, mu

C material properties
      Y = props(1)      ! Young's modulus
      n = props(2)      ! Poisson's ratio

C Lame's parameters
      lambda = Y*n/((1.0d0+n)*(1.0d0-2.0d0*n))
      mu = Y/(2.0d0*(1.0d0+n))

C Stiffness matrix
      do i = 1, ntens
        do j = 1, ntens
          ddsdde(i,j) = 0.0d0
        end do
      end do
      do i = 1, ndi
        do j = 1, ndi
          ddsdde(i, j) = lambda
        end do 
        ddsdde(i,i) = lambda + 2.0d0*mu
      end do 

C Shear contribution
      do i = ndi+1, ntens
        ddsdde(i,i) = mu
      end do 

C Stress increment evaluation
      do i = 1, ntens
        do j = 1, ntens
          stress(i) = stress(i) + ddsdde(i,j) * dstran(j)
        end do 
      end do 
C
      return
      end