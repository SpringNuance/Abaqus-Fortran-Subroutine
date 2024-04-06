   include 'vaba_param.inc'
   parameter( nrData=6 )  ! maximum number of components of variables to be requested (6 for symmetric tensor)
   character*3 cData(maxblk*nrData)  ! parameter maxblk is the maximum value of nblock
   dimension rData(maxblk*nrData), jData(maxblk*nrData)
   ...
   call vgetvrm( 'VAR', rData, jData, cData, jStatus )