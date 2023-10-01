! -----------------------------------------------------------!
! Subroutine for Abaqus/Explicit for isotropic elasticity    !
! and isotropic plasticity with pressure and lode            !
! dependence, according to plasticity model by BAI, Y. and   !
! WIERZBICKI, T. (published in "A new model of metal         !
! plasticity and fracture with pressure and Lode             !
! dependence", International Journal of Plasticity 24 (2008) !
! pages 1071-1096).                                          !
!						             !
! - Damage initiation according to 3D fracture locus         !
!   (equivalent plastic strain to failure initiation being   !
!   function of both hydrostatic pressure and lode angle)    !
! - Damage evolution based on energy dissipated during       !
!   damage process                                           !
!							     !
! - Temperature softening				     !
! - Strain rate hardening				     !
!                                                            !
! Originally based on the von Mises' J2 plasticity theory.   !
! Elastic predictor, radial corrector algorithm used. With   !
! explicit forward Euler scheme for integration of flow rule.!
!                                                            !
! Not suitable for solving problem settings including plane  !
! stresses. Not suitable for 2D models.                      !
! -----------------------------------------------------------!
! Solution dependent variables (SDVs):			     !
!							     !
! SDV1: 	Equivalent plastic strain		     !
! SDV2: 	Damage variable									 
! SDV3: 	Yield stress at damage initiation				 
! SDV4: 	Flag (=0 element not damaged,					 
!			      =1 element experienced brittle fracture,	 
!			      =2 element experienced ductile damage,     
!                 =3 element experienced ductile fracture)	 
! SDV5-10:	Total strain tensor								 
! SDV11-16:	Plastic strain tensor							 
! SDV17:	Equivalent plastic strain increment				 
! SDV18:	Temperature softening correction function		 
! SDV19:	Mutliplicative strain rate hardening correction	 
!			function		 								 !
! SDV20:	Equivalent plastic strain rate					 
! SDV21:	Equivalent plastic strain rate (computed over	 
!           a user-defined no. of increments)				 
! SDV22:	Increment counter (for eq. plas. str. rate       
!           computation)									 
! SDV23:	Sum of equivalent plastic strain increments      
!           before the computation of eq. plas. str. rate    
! SDV24:	Sum of time increments before computation of eq. 
!           plas. str. rate									 
! SDV25:	Additive strain rate hardening correction		 
!			function										 !
! SDV26:	Maximum principal stress						 
! SDV27:	Equivalent plastic strain (averaged over a		 
!			user-defined no. of increments)					 
! SDV28:	Element deletion flag							 
!------------------------------------------------------------!
! SDV and parameters modified by B. Wu, J. He and F. Shen    !
!                                                            !
! SDV29:	Eta         								     
! SDV30:	Thetabar										 !
! SDV31:	Ductile damage initiation indicator				 
! SDV32:	Ductile failure indicator						 
! SDV33:	Damage initiation strain at the current step	 
! SDV34:	Dcrit at the current step						 
! SDV35:	Cleavage initiation strain at the current step	 
! SDV36:	Cleavage failure indicator based on strain		 
! SDV37:	Damage fracture strain at the current step	     
! SDV38:	Equivalent plastic strain at damage initiation	 
! SDV39:	Temperatrue	 								     !
! -----------------------------------------------------------!
! Additional parameters props(28)-props(36)                  !
!                                                            !
! Additional parameter for cut-off value is defined	     !
! Additional parameters Ddf1-Ddf4 are defined for Ddf-locus  !
! Additional parameters Cc1-Cc4 are defined for Cleavage     !
! -----------------------------------------------------------!
! Modifications by B. Wu, J. He and F. Shen                  !
!                                                            !
! Damage indicator: Damage initiates when SDV31>1            !
! Failure indicator: Failure happens when SDV32>1            !
! Cleavage failure indicator: Cleavage happens when SDV36>1  !
! Non-proportional Loading paths with varying stress states  !
! are considered are considered with weighting scheme        !
! Cut-off value concept is considered                        !
! -----------------------------------------------------------!

SUBROUTINE VUMAT( &
! Read only
     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
     stepTime, totalTime, dt, cmname, coordMp, charLength, &
     props, density, strainInc, relSpinInc, &
     tempOld, stretchOld, defgradOld, fieldOld, &
     stressOld, stateOld, enerInternOld, enerInelasOld, &
     tempNew, stretchNew, defgradNew, fieldNew, &
! Write only
     stressNew, stateNew, enerInternNew, enerInelasNew)
     
     include 'vaba_param.inc'

        REAL :: props(nprops)
        REAL :: density(nblock)
        REAL :: coordMp(nblock,:)
        REAL :: charLength(nblock)
        REAL :: strainInc(nblock,ndir+nshr)
        REAL :: relSpinInc(:)
        REAL :: tempOld(nblock)
        REAL :: stretchOld(:)
        REAL :: defgradOld(:)
        REAL :: fieldOld(:)
        REAL :: stressOld(nblock,ndir+nshr)
        REAL :: stateOld(nblock,nstatev)
        REAL :: enerInternOld(nblock)
        REAL :: enerInelasOld(nblock)
        REAL :: tempNew(nblock)
        REAL :: stretchNew(:)
        REAL :: defgradNew(:)
        REAL :: fieldNew(:)
        REAL :: stressNew(nblock,ndir+nshr)
        REAL :: stateNew(nblock,nstatev)
        REAL :: enerInternNew(nblock)
        REAL :: enerInelasNew(nblock)


        character*80 cmname

! Defining numerical constants
        parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     1  threeHalves = 1.5, thousand = 1000., tol=1.D-6 )

! Defining material properties constants
        e      = props(1)
        xnu    = props(2)
        ceta   = props(3)
        eta0   = props(4)
        cthetas= props(5)
        cthetat= props(6)
        cthetac= props(7)
        om     = props(8)
! Ductile damage initiation locus function parameters
        Ddi1   = props(9)
        Ddi2   = props(10)
        Ddi3   = props(11)
        Ddi4   = props(12)
        Ddi5   = props(13)
        Ddi6   = props(14)
! Fracture energy dissipation (energy dissipated per 
! element unit area after ductile damage dissipation)
        gf     = props(15)
! Maximum value for ductile damage variable
!       Dcrit   = props(16)
! 5 parameters for temperature softening correction
        cT1    = props(17)
        cT2    = props(18)
        cT3    = props(19)
        eta2   = props(20)
        cp     = props(21)
        t0     = props(22)
! 2 parameters for strain rate hardening correction
        cE1    = props(23)
        cE2    = props(24)
        cE3    = props(25)
! Number of increments to update plastic strain rate
        strrInc= props(26)
! Brittle damage initiation stress value
        sigdmg = props(27)
! Cut-off value
        Cut    = props(28)
! Ductile frature strain locus function parameters
        Ddf1   = props(29)
        Ddf2   = props(30)
        Ddf3   = props(31)
        Ddf4   = props(32)
! Cleavage initiation locus
        Cc1    = props(33)
        Cc2    = props(34)
        Cc3    = props(35)
        Cc4    = props(36)		
		
        nvalue = (nprops-40)/2
!
! Checking for valid entries
        if (om.lt.zero) then
          write(6,5)
 5        format(//,30X,'***ERROR - m MUST BE A NON-NEGATIVE INTEGER')
        endif
!
! Defining constants
        twomu  = e / ( one + xnu )
        thremu = threeHalves * twomu
        sixmu  = three * twomu
        alamda = twomu * ( e - twomu ) / ( sixmu - two * e )
        akappa = e * twomu * half / (three * (threeHalves * twomu - e))
        term   = one / ( twomu * ( one + hard/thremu ) )
        con1   = sqrt( twoThirds )
        pi     = 4. * atan(one)
        con2   = (sqrt(three)/two) / (one - sqrt(three)/two)
!
! Computation per material point starts here
      do 100 i = 1,nblock
!
! If brittle fracture has previously occured, stress is zeroed
        if (stateOld(i,4).gt.half .and. stateOld(i,4).lt.(one+half))
     1    then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,2) = one
		  stateNew(i,4) = one
		  stateNew(i,28) = zero
          goto 100
        endif
! If ductile fracture has previously occured, stress is zeroed
        if (stateOld(i,4).gt.(two+half) .and. stateOld(i,4).lt.(three+half))
     1    then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
		  stateNew(i,2) = one
          stateNew(i,4) = three
		  stateNew(i,32) = one
		  stateNew(i,28) = zero
          goto 100
        endif
! Updating total strain (states 5-10)
        stateNew(i,5)  = stateOld(i,5)  + strainInc(i,1)
        stateNew(i,6)  = stateOld(i,6)  + strainInc(i,2)
        stateNew(i,7)  = stateOld(i,7)  + strainInc(i,3)
        stateNew(i,8)  = stateOld(i,8)  + strainInc(i,4)
        stateNew(i,9)  = stateOld(i,9)  + strainInc(i,5)
        stateNew(i,10) = stateOld(i,10) + strainInc(i,6)

! Trace of total strain tensor
        epsilontrace = stateNew(i,5) + stateNew(i,6) + stateNew(i,7)

! Trace of strain increment tensor
        trace = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)

! Calculating softening correction term due to temperature rise
        if (tempOld(i).le.tol .and. tempOld(i).ge.-tol) then
          facT = cT1*exp(zero-cT2*t0)+cT3
          stateNew(i,18)= facT
		  stateNew(i,39)= t0
        else
          facT = cT1*exp(zero-cT2*tempOld(i))+cT3
          stateNew(i,18)= facT
		  stateNew(i,39)= tempOld(i)
        endif   

! Calculating hardening correction terms due to straining rate
        facStrrateM = cE1*log(stateOld(i,21)) + cE2
        facStrrateA = cE3*stateOld(i,21)

! Preventing negative values for straining rate
        if (facStrrateM.lt. one) facStrrateM = one
        stateNew(i,19)= facStrrateM

! Restoring previous softening correction term due to ductile damage
        facDctlDmgDev = one - stateOld(i,2)
        facDctlDmgVol = facDctlDmgDev
        if (stateOld(i,29).lt.zero) then
          facDctlDmgVol = one
        endif

! Trial stress
        sig1 = stressOld(i,1) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,1)
        sig2 = stressOld(i,2) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,2)
        sig3 = stressOld(i,3) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,3)
        sig4 = stressOld(i,4) + facDctlDmgDev*twomu*strainInc(i,4)
        sig5 = stressOld(i,5) + facDctlDmgDev*twomu*strainInc(i,5)
        sig6 = stressOld(i,6) + facDctlDmgDev*twomu*strainInc(i,6)

! Following equation numbering of publication mentioned in code title text
! Equ. (4) - Calculating the deviatoric part of trial stress
        sigmean = third * ( sig1 + sig2 + sig3 )
        ds1 = sig1 - sigmean
        ds2 = sig2 - sigmean
        ds3 = sig3 - sigmean

! Calculating the magnitude of the deviatoric trial stress tensor
        dsmag = sqrt( ds1**2 + ds2**2 + ds3**2 + two*sig4**2 + two
     1   *sig5
     1   **2 + two*sig6**2 )

! Preventing a divide by zero when dsmag is zero. When dsmag is zero, computation is still in the elastic zone
        if (dsmag.lt.tol .and. dsmag.ge.zero) dsmag = tol
        if (dsmag.gt.(zero-tol) .and. dsmag.le.zero) dsmag = zero-tol

! Following equation numbering of publication mentioned in code title text
! Eq. (1) - Calculating the 1st invariant of the stress tensor
        p   = zero - sigmean

! Eq. (2) - Calculating the 2nd invariant of the stress tensor
        q   = dsmag / con1

! Eq. (3) - Calculating the 3rd invariant of the stress tensor
        r   = ( three/twoThirds * (ds1*(ds1**2 + sig4**2 + sig6**2)
     1  + 2.
     1  *sig4*(ds1*sig4 + ds2*sig4 + sig6*sig5) + 2.*sig6*(ds1*sig6
     1  + sig4
     1  *sig5 + ds3*sig6) + ds2*(sig4**2 + ds2**2 + sig5**2) + 2.*sig5
     1  *(sig4*sig6 + ds2*sig5 + ds3*sig5) + ds3*(sig6**2 + sig5**2
     1  + ds3**2)) )**third

! Eq. (5) - Calculating the dimensionless hydrostatic pressure eta
        eta = sigmean / q
        stateNew(i,29)= eta

! Eq. (6) - Calculating the Lode angle theta
        cosine=(r/q)**3
!       Assuring that -1<cosine<1
        if (cosine.gt.one) cosine=one
        if (cosine.lt.(zero-one)) cosine=zero-one
        theta = third*acos(cosine)

! Eq. (7) - Calculating the normalized Lode angle thetabar
        thetabar = one - 6.*theta/pi
        stateNew(i,30)= thetabar

! Eq. (12) - Calculating the Lode angle dependant parameter gamma
        gamma = con2 * (one/cos(theta-pi/6.) - one)

! Eq. (13) - Determining cthetaax, either tension case or compression case
        if( thetabar .ge. zero ) then
          cthetaax = cthetat
        else
          cthetaax = cthetac
        endif

! Fetching yield stress from flow curve
        call ahard(sigmayield,hard,stateOld(i,1),props(41),nvalue)

! Eq. (11) - Calculating radius of elastic zone
        radius = (sigmayield * facStrrateM + facStrrateA) * ( one
     1  - ceta
     1	*(eta-eta0) ) * (cthetas + ((cthetaax-cthetas)* (gamma- ((gamma
     1  **(om+one))/(om+one)) ) ) ) * facT * facDctlDmgVol

! Eq. (15) - Checking if yielding. The yielding factor facyld is set to be zero
! in the elastic zone, one if yielding
        facyld = zero
        if( q - radius .ge. zero ) facyld = one

! Eq. (20) - Calculating the tensor of the normal direction of plastic flow
!       Avoiding a divide by zero by avoiding that theta = zero
        if (theta.eq.zero) theta = theta + tol

        fac1 = con1 * three/(two*q)
        fac2 = facT * facDctlDmgVol * con1 * (sigmayield * facStrrateM
     1  + facStrrateA)*ceta*(cthetas+((cthetaax-cthetas)*(gamma
     1  -( (gamma**(om+one)) / (om+one) )))) * three*eta/(two*(q**2))
        fac3 = facT * facDctlDmgVol * con1 * (sigmayield * facStrrateM
     1  + facStrrateA) * (one-ceta*(eta-eta0)) * (cthetaax-cthetas)
     1  * (one -(gamma**om)) * (three*sqrt(three)/(two-sqrt(three)))
     1  * (tan(theta-pi/6.)/cos(theta-pi/6.)) *(one/(q*sin(three
     1  *theta)))

        on1 = fac1*ds1
     1  - fac2 * ds1
     1  - fac3 * (one/three+cos(three*theta)*ds1/(two*q)-three*((ds1
     1  **2)+(sig4**2)+(sig6**2))/(two*(q**2)))

        on2 = fac1*ds2
     1  - fac2 * ds2
     1  - fac3 * (one/three+cos(three*theta)*ds2/(two*q)-three*((sig4
     1  **2)+(ds2**2)+(sig5**2))/(two*(q**2)))

        on3 = fac1*ds3
     1  - fac2 * ds3
     1  - fac3 * (one/three+cos(three*theta)*ds3/(two*q)-three*((sig6
     1  **2)+(sig5**2)+(ds3**2))/(two*(q**2)))

        on4 = fac1*sig4
     1  - fac2 * sig4
     1  - fac3 * (cos(three*theta)*sig4/(two*q)-three*((ds1*sig4)
     1  +(ds2*sig4)+(sig6*sig5))/(two*(q**2)))

        on5 = fac1*sig5
     1  - fac2 * sig5
     1  - fac3 * (cos(three*theta)*sig5/(two*q)-three*((sig4*sig6)
     1  +(ds2*sig5)+(ds3*sig5))/(two*(q**2)))

        on6 = fac1*sig6
     1  - fac2 * sig6
     1  - fac3 * (cos(three*theta)*sig6/(two*q)-three*((ds1*sig6)
     1  +(sig4*sig5)+(ds3*sig6))/(two*(q**2)))

! Calculating trial yield function
        phitrial   = facyld * con1 * (q - radius) / facDctlDmgVol

! Calculating equivalent plastic strain
        deqps  = con1 * phitrial / (twomu + twoThirds * hard)
        if (deqps.lt.zero) then
			stateNew(i,1) = stateOld(i,1)
			goto 100
        endif
		stateNew(i,17) = deqps

! Updating equivalent plastic strain
        stateNew(i,1) = stateOld(i,1) + deqps

! Updating plastic strain (states 11-16)
        stateNew(i,11) = stateOld(i,11) + deqps * on1 / con1
        stateNew(i,12) = stateOld(i,12) + deqps * on2 / con1
        stateNew(i,13) = stateOld(i,13) + deqps * on3 / con1
        stateNew(i,14) = stateOld(i,14) + deqps * on4 / con1
        stateNew(i,15) = stateOld(i,15) + deqps * on5 / con1
        stateNew(i,16) = stateOld(i,16) + deqps * on6 / con1

! Updating equivalent plastic strain rate
        stateNew(i,22) = stateOld(i,22) + one
        stateNew(i,23) = stateOld(i,23) + deqps
        stateNew(i,24) = stateOld(i,24) + dt
        if (strrInc.gt.half .and. strrInc.lt.(one+half)) then
			stateNew(i,21) = deqps/dt
			stateNew(i,27) = deqps
        else
          stateNew(i,21) = (stateNew(i,23)/stateNew(i,24)
     1                   +  stateOld(i,21))/two
          stateNew(i,27) = (stateNew(i,23)/stateNew(i,22)
     1                   +  stateOld(i,27))/two
          if (stateNew(i,22).gt.strrInc-half .and. stateNew(i,22)
     1    .lt.strrInc+half) then
            stateNew(i,22) = zero
            stateNew(i,23) = zero
            stateNew(i,24) = zero
          endif
        endif
        stateNew(i,20) = deqps / dt

! Updating temperature
        if (tempOld(i).le.tol .and. tempOld(i).ge.-tol) then
          tempNew(i) = t0
     1    +(eta2/(density(i)*cp))*radius*deqps/facDctlDmgDev
        else
          tempNew(i) = tempOld(i)
     1    +(eta2/(density(i)*cp))*radius*deqps/facDctlDmgDev
        endif

! Calculating equivalent plastic strain to damage initiation
        epsilondi = (half * ( Ddi1*exp(zero-Ddi2*eta) + Ddi5*exp(zero-Ddi6*eta))
     1                     - Ddi3*exp(zero-Ddi4*eta) )
     1                     * thetabar**2
     1           +  half * ( Ddi1*exp(zero-Ddi2*eta) - Ddi5*exp(zero-Ddi6*eta))
     1                     * thetabar
     1           +  Ddi3*exp(zero-Ddi4*eta)
! Jimmy delete          stateNew(i,33)= epsilondi
!  update damage initiation indicator, damage initition is suppressed when 
!  the stress triaxiality is less than the cut-off value
        if (eta.lt. Cut .and. stateNew(i,31).lt.1.) then
          stateNew(i,31) = stateOld(i,31)
        else
        stateNew(i,31) = stateOld(i,31)+deqps/epsilondi
        endif
        if (stateNew(i,31).ge.1.) then
			stateNew(i,31) = 1.0
        endif

! Calculating equivalent plastic strain to cleavage failure
       epsilonc=  ( Cc1*exp(zero-Cc2*eta) - 
     1                    Cc3*exp(zero-Cc4*eta))* thetabar**2
     1           +  Cc3*exp(zero-Cc4*eta)
       stateNew(i,35) = epsilonc
	   stateNew(i,36) = stateOld(i,36)+ deqps/epsilonc	
        if (stateNew(i,36).ge.1.) then
			stateNew(i,36) = 1.0
        endif
	   
! Calculating equivalent plastic strain to ductile failure
       epsilondf=  (Ddf1*exp(zero-Ddf2*eta) - 
     1                    Ddf3*exp(zero-Ddf4*eta))* thetabar**2
     1           +  Ddf3*exp(zero-Ddf4*eta)
! Jimmy change
       delta_df=deqps/(epsilondf-epsilondi)	 
! Jimmy delete       stateNew(i,37) = epsilondf   

! SDV 27: deqps
! SDV 31: Idd
! SDV 32: Idf
! SDV 33: PEEQ@ddi
! SDV 34: Dcr
! SDV 35: 
! SDV 36: 
! SDV 37: PEEQ@df
! SDV 38: PEEQ@ddi the same as 33
! SDV 39: 
! 
! Ductile damage

        if (stateNew(i,31).ge.1.0 .and. stateNew(i,1).gt.tol) then
          if (stateOld(i,3).le.tol .and. stateOld(i,3).ge.-tol) then
!           Registering the yield stress at the onset of damage for the first time
            stateNew(i,3) = radius
!           Jimmy change
			stateNew(i,33) = epsilondi
			stateNew(i,37) = epsilondf
!
			stateNew(i,38) = stateNew(i,1)
          else
!           Keeping the stress and PEEQ at the onset of damage for next step
            stateNew(i,3) = stateOld(i,3)
!           Jimmy change
			stateNew(i,33) = stateOld(i,33)
			stateNew(i,37) = stateOld(i,37)
!
			stateNew(i,38) = stateOld(i,38)
          endif
!         Updating failure indicator and damage fraction variable
          Dcrit=(stateNew(i,3))*(stateNew(i,37)-stateNew(i,33))/(gf)
          stateNew(i,34)= Dcrit
          if (eta.gt. Cut) then
            if (stateOld(i,32).lt.one) then
!               stateNew(i,32) = stateOld(i,32) + deqps/epsilondf
               stateNew(i,32) = stateOld(i,32) + delta_df
!			    stateNew(i,2) = stateOld(i,2) + Dcrit*(deqps/epsilondf)
			   stateNew(i,2) = Dcrit*stateNew(i,32)*zero
            else
             stateNew(i,2) = one
			 stateNew(i,32) = one
			 stateNew(i,4) = three
			 stateNew(i,28) = zero
            endif
          endif
!       Elements under hydrostatic compression don't experience spherical damage
          if (eta.lt. Cut .and. stateNew(i,32).lt.one) then
            stateNew(i,2) = stateOld(i,2)
            stateNew(i,32) = stateOld(i,32)
          endif
!         Registering color code for elements
          if(stateNew(i,31).ge.1.0 .and. stateNew(i,32).lt.one) then
		    stateNew(i,4) = two
          else
		    stateNew(i,4) = three
          endif
        else
!         In case no damage is reached, yield stress at damage onset
!         and damage variable are kept for next step
          stateNew(i,3) = stateOld(i,3)
		  stateNew(i,38) = stateOld(i,38)
          stateNew(i,2) = stateOld(i,2)
		  stateNew(i,32) = stateOld(i,32)
        endif

!       Previous ductile damage
        if (stateOld(i,4).gt.(one+half)
     1  .and. stateOld(i,4).lt.(two+half)) then
          stateNew(i,4) = two
        endif
        if (stateOld(i,4).gt.(two+half)
     1  .and. stateOld(i,4).lt.(three+half)) then
          stateNew(i,4) = three
        endif
!       Jimmy change avoid D larger than one situation
        if (stateNew(i,2).lt.one) then		
!       	Calculating new softening correction terms due to ductile damage
			stateNew(i,28)= one
			if (stateNew(i,32).ge.one) then
			  facDctlDmgDev = zero
			  facDctlDmgVol = zero
			  stateNew(i,32) = one
			  stateNew(i,2) = one
			  stateNew(i,4) = three
			  stateNew(i,28)= zero
			  goto 100
			else
			  facDctlDmgDev = one - stateNew(i,2)
			  facDctlDmgVol = facDctlDmgDev
			endif
		else
			facDctlDmgDev = zero
			facDctlDmgVol = zero
			stateNew(i,32) = one
			stateNew(i,2) = one
			stateNew(i,4) = three
			stateNew(i,28)= zero
		endif
! 		Jimmy change: avoid extensive distortion 
!		dist_ratio = stataeNew(i,27)/stateOld(i,27) 
!		if (dist_ratio.ge.(two+three)) then 
!			facDctlDmgDev = zero
!			facDctlDmgVol = zero
!			stateNew(i,32) = one
!			stateNew(i,2) = one
!			stateNew(i,4) = three
!			stateNew(i,28)= zero
!		endif
!       Elements under hydrostatic compression don't experience spherical damage
        fac1 = facDctlDmgVol * akappa * epsilontrace
        fac2 = facDctlDmgDev * twomu * deqps / con1
        sig1 = facDctlDmgDev * twomu * (stateNew(i,5)
     1       - epsilontrace/three - stateOld(i,11))
        sig2 = facDctlDmgDev * twomu * (stateNew(i,6)
     1       - epsilontrace/three -stateOld(i,12))
        sig3 = facDctlDmgDev * twomu * (stateNew(i,7)
     1       - epsilontrace/three -stateOld(i,13))
        sig4 = facDctlDmgDev * twomu * (stateNew(i,8) -stateOld(i,14))
        sig5 = facDctlDmgDev * twomu * (stateNew(i,9) -stateOld(i,15))
        sig6 = facDctlDmgDev * twomu * (stateNew(i,10)-stateOld(i,16))
! Update stress		
        stressNew(i,1) = sig1 + fac1 - fac2 * on1
        stressNew(i,2) = sig2 + fac1 - fac2 * on2
        stressNew(i,3) = sig3 + fac1 - fac2 * on3
        stressNew(i,4) = sig4        - fac2 * on4
        stressNew(i,5) = sig5        - fac2 * on5
        stressNew(i,6) = sig6        - fac2 * on6
!       Calculating invariants of stress tensor
        sig1 = stressNew(i,1)
        sig2 = stressNew(i,2)
        sig3 = stressNew(i,3)
        sig4 = stressNew(i,4)
        sig5 = stressNew(i,5)
        sig6 = stressNew(i,6)
        SI1 =   sig1 + sig2 + sig3
        SI2 =   sig1*sig2-sig4*sig4
     1        + sig1*sig3-sig6*sig6
     1        + sig2*sig3-sig5*sig5
        SI3 =   sig1*(sig2*sig3-sig5*sig5)
     1        - sig4*(sig4*sig3-sig5*sig6)
     1        + sig6*(sig4*sig5-sig2*sig6)
!       Preparing subvalues for calculating the principal stresses values
        cosine2 = (two*SI1*SI1*SI1-three*three*SI1*SI2+three*three
     1            *three*SI3)
     1            /(two*(SI1*SI1-three*SI2)**(threehalves))
!       Assuring that -1<cosine2<1
        if (cosine2.gt.one) cosine2=one
        if (cosine2.lt.(zero-one)) cosine2=zero-one
        alpha2 = acos(cosine2)
!       Calculating the principal stress values
        SP1 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three)
        SP2 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three+twoThirds*pi)
        SP3 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(twoThirds*pi-alpha2/three)
!       Fetching the highest of the principal stress values
        sigmamax = max(abs(SP1),abs(SP2),abs(SP3))
        stateNew(i,26) = sigmamax

!       Brittle damage
      if (stateNew(i,36).ge.one .and. stateNew(i,1).gt.tol) then
        if (sigmamax.ge.sigdmg .and. stressNew(i,1).gt.zero) then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,4)  = one
		  stateNew(i,2)  = one
          stateNew(i,28) = zero
           write(*,*) "Coordinate:    ", coordMp(i,1), coordMp(i,2), coordMp(i,3)
           write(*,*) "stateNew(i,1): ", stateNew(i,1), epsilonc
           write(*,*) "eta thetabar: ", eta, thetabar
           write(*,*) "maximum principal stress:", sigmamax
        endif
      endif		
		
  100 continue
!
      return
      end

      subroutine ahard(sigmayield,hard,eqplas,table,nvalue)

      include 'vaba_param.inc'
      dimension table(2,nvalue)

!     Set yield stress to last value of table, hardening to zero
      sigmayield=table(1,nvalue)
      hard=0.0

!     If more than one entry, search table
      if(nvalue.gt.1) then
        do 10 k1=1,nvalue-1
          eqpl1=table(2,k1+1)
          if(eqplas.lt.eqpl1) then
            eqpl0=table(2,k1)
            if(eqpl1.le.eqpl0) then
              write(6,7)
 7            format(//,30X,'***ERROR - PLASTIC STRAIN MUST BE ',
     1               'ENTERED IN ASCENDING ORDER,')

!             Subroutine XIT terminates execution and closes all files
!              call XIT
            endif
            deqpl=eqpl1-eqpl0
            sigmayield0=table(1,k1)
            sigmayield1=table(1,k1+1)
            dsigmayield=sigmayield1-sigmayield0
            hard=dsigmayield/deqpl
            sigmayield=sigmayield0+(eqplas-eqpl0)*hard
            goto 20
          endif
 10     continue
 20     continue
        if(eqplas.gt.table(2,nvalue)) then
          hard=(table(1,nvalue)-table(1,nvalue-1))
     1        /(table(2,nvalue)-table(2,nvalue-1))
          sigmayield=table(1,nvalue)+(eqplas-table(2,nvalue))*hard
        endif
      endif
      return

! Iteration ends here
      end