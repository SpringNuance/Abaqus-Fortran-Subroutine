c ======================================================================
c User Subroutine UMAT for Abaqus linear elastic material
c By Irfan Habeeb CN (PhD, Technion - IIT)
c ======================================================================
      subroutine vumat(
C Read only -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew)
C
      include 'vaba_param.inc'
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

      character*80 cmname
      integer k, k1
C inputs
      real*8 e, nu, mu, alamda, Scrit, mises, S, Strace, damage

      e = props(1)                  ! Young's modulus
      nu = props(2)                 ! poisson's ratio
      rho = props(3)                ! density
      Scrit = props(4)              ! critical stress

c lame's parameters
      mu = e/(2.d0*(1.d0 + nu))
      alamda = e*nu/((1.d0 + nu) * (1.d0 - 2.d0*nu))

c stress increment evaluation for each element
      do 10 k = 1, nblock
        trInc = sum(strainInc(k, 1:3))        ! strain trace
        do k1 = 1, ndir
          stressNew(k, k1) = stressOld(k, k1) + alamda*trInc 
     1      + 2.d0*mu*strainInc(k, k1)
        end do
c shear stress
        do k1 = 1, nshr
          stressNew(k, ndir+k1) = stressOld(k, ndir+k1) + 
     1      2.d0*mu*strainInc(k, ndir+k1)
        end do

c Von mises stress
        Strace = sum(stressOld(k, 1:ndir))
        S = 0.d0
        do k1 = 1, ndir
          S = S + (stressOld(k, k1) - Strace)**2.d0
        end do 
        do k1 = ndir+1, ndir+nshr
          S = S + 2.d0 * stressOld(k, k1)**2.d0
        end do 
        mises = sqrt(3.d0/2.d0 * S)

c defining the damage criteria, if the mises exceeds the critical stress
        damage = 1.d0
        if ((mises .gt. Scrit) .and. (stateOld(k, 2).gt.0.d0)) then
          damage = 0.d0
        end if
        
c updating the state variables
        stateNew(k, 1) = damage         ! element damage
c hydrostatic strain, to constrain damage only under tension
        stateNew(k, 2) = stateOld(k, 2) + sum(strainInc(k, 1:ndir))/3.d0

c updating the internal energy
      enerInternNew(k) = enerInternOld(k) + 
     1 dot_product( 0.5d0*(stressOld(k, 1:ndir+nshr) +
     2 stressNew(k, 1:ndir+nshr)), strainInc(k, 1:ndir+nshr)) / rho

 10   continue
      return
      end