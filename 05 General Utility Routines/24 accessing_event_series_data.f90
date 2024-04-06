include 'aba_evs_param.inc'

character*80 evsName, cProps(n_evsC_size)
dimension iProps(n_evsI_size), rProps(n_evsR_size)
...
call getEventSeriesSliceProperties(evsName, iProps, rProps, cProps)  


include 'aba_evs_param.inc'

character*80 evsName
dimension  evsTime(maxEvents), evsCoord(3,maxEvents), evsField(numFields,maxEvents)
...
call getEventSeriesSliceLG(evsName, numEvsPoints, evsTime, evsCoord, evsField) 


include 'aba_evs_param.inc'

character*80 evsName
dimension  xCenter(3), evsPathTime(2,maxEvents), evsPathCoord(3,2,maxEvents),
 evsPathField(nFieldSize,2,maxEvents)
...
call getEventSeriesSliceLGLocationPath(evsName, xCenter, radius, 
nFieldSize, numSegments, evsPathTime, evsPathCoord, evsPathField) 