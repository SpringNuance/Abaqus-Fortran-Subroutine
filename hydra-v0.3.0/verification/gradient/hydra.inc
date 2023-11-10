      PARAMETER (nelhmax=16384)
      PARAMETER (nndhmax=16384)

C     CPS4R
      PARAMETER (ND=2)
      PARAMETER (NELHINC=4)
      PARAMETER (NELHGPT=1)

      COMMON /kELHI/ nelh, kelhlab(nelhmax), kelhinc(NELHINC,nelhmax)
      COMMON /kELHD/ delhgradph(ND,NELHGPT,nelhmax),
     +               delhjac(ND,NELHINC,NELHGPT,nelhmax)
