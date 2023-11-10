      PARAMETER (nelhmax=16384)
      PARAMETER (nndhmax=16384)
      PARAMETER (nelcohmax=1000)

C     CPx4R{T}
      PARAMETER (ND=2)
      PARAMETER (NELHINC=4)
      PARAMETER (NELHGPT=1)

      COMMON /kELHI/ nelh, kelhlab(nelhmax), kelhinc(NELHINC,nelhmax)
      COMMON /kELHD/ delhgradph(ND,NELHGPT,nelhmax), 
     +               delhjac(ND,NELHINC,NELHGPT,nelhmax),
     +               delhpeeq0(NELHGPT,nelhmax),
     +               delhpeeqC(NELHGPT,nelhmax),
     +               delhpeeqCtime(NELHGPT,nelhmax)


C     COH2D4 state variables indexes
      PARAMETER (N_DAM=1)
      PARAMETER (N_coh_SDV=1)

C     AISI material user defined fields
      PARAMETER (NFIELD_AISI=1)
      PARAMETER (IPHI_AISI=1)
C     AISI state variables indexes
      PARAMETER (I_aNt=1, I_Cl=2, I_PEEQ=3, I_Ct=4, I_C=5, I_theta=6,
     +           I_dec=7, I_phi=8)
      PARAMETER (NSTATV_AISI=8)


      COMMON /kELCOI/ nelcoh, kelmap(2,nelcohmax)
      COMMON /kELCOD/ dcohdam(nelcohmax)
