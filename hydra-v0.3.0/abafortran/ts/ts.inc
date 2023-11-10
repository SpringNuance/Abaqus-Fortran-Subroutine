c vi:syntax=fortran

c component ordering in arrays returned by GETVRM
c see AUSRG 2.1.6 'Ordering of returned components'
      PARAMETER (I_11=1, I_22=2, I_33=3)

c variable ordering in array dprop
      PARAMETER (i_LE11=1)

c static dimensioning of common blocks
      PARAMETER (nelcohmax=2000)
      PARAMETER (nmap=2, nprop=1)

      COMMON /kCOST/nelmap
      COMMON /kMAP/kelmap(nmap,nelcohmax)
      COMMON /kPROP/dprop(nprop,nelcohmax)

c kelmap(1:2,1:nelmap) = ((CPE4R, COH) ... )
c dprop(1:1,1:nelmap) = ((LE11) ... )
