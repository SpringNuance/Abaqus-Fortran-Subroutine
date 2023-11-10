      PARAMETER (nelcohmax=2000)
      PARAMETER (N_DEC=1,N_NNC=2)

      COMMON /kCOST/nelmap
      COMMON /kMAP/kelmap(2,nelcohmax)
      COMMON /kPROP/dprop(1,nelcohmax)

C   kelmap(1:2,1:nelmap) = ((CPE4R, COH) ... )
