character*80 tcName 
...
call setTableCollection(tcName, jError)


include 'aba_tcs_param.inc'

character*80 parameterTableLabel, cParams(maxParams)
dimension iParamsDataType(maxParams), iParams(maxParams), rParams(maxParams)
...
call getParameterTable(parameterTableLabel, numParams, iParamsDataType, iParams, rParams, cParams, jError)



include 'aba_tcs_param.inc'

character*80 parameterTableLabel, cParams(maxParams)
dimension iParamsDataType(maxParams), iParams(maxParams), rParams(maxParams)
...
call getParameterTableRow(parameterTableLabel, jRow, numParams, iParamsDataType, iParams, 
rParams, cParams, jError)


character*80 propertyTableLabel
dimension field(*)
dimension rIndepVars(maxIndepVars), rProps(maxProps), dPropDVar(maxProps*maxVars) 
...
call getPropertyTable(propertyTableLabel, rIndepVars, temp, field, 
                                      numProps, rProps, dPropDVar, numDervs, jError)


character*80 propertyTableLabel
dimension temp(nBlock), field(nBlock,*)
dimension rIndepVars(maxBlock*maxIndepVars), rProps(maxBlock*maxProps), 
          dPropDIVar(maxBlock*maxProps*maxIndepVars), dPropdTemp(maxBlock*maxProps)
...
call vGetPropertyTable(propertyTableLabel, nBlock, rIndepVars, temp, field, 
                      numProps, rProps, dPropDIVar, dPropDTemp, jDervs, jError)


character*80 propertyTableLabel
integer*8 mIndepVars(maxIndepVars)
dimension temp(nBlock), field(nBlock,*)
dimension rProps(maxBlock*maxProps), dPropDIVar(maxBlock*maxProps*maxIndepVars), dPropdTemp(maxBlock*maxProps)
...
call SMASetPointer(mIndVarPtrs(1), eqps)
call SMASetPointer(mIndVarPtrs(2), eqpsRate)
call vGetPropertyTablePtr(propertyTableLabel, nBlock, mIndepVars, temp, field, 
                      numProps, rProps, dPropDIVar, dPropDTemp, jDervs, jError)


include 'aba_tcs_param.inc'

dimension jSize(n_tcsI_TCDB_size)
...
call queryTableCollectionSize(jSize, jError)



character*80 tcNames(maxTCs)
...
call queryTableCollectionNames(tcNames, numTCs)


include 'aba_tcs_param.inc'

character*80 cTableColl(n_tcsC_TC_size)
dimension jTableColl(n_tcsI_TC_size)
...
call queryTableCollection(jTableColl, cTableColl, jError)


character*80 parameterTableLabel
...
queryParameterTable(parameterTableLabel, numParams, numRows, jError)


include 'aba_tcs_param.inc'

character*80 propertyTableLabel, cPropTable(n_tcsC_PRT_size)
dimension jPropTable(n_tcsI_PRT_size), rPropTable(n_tcsR_PRT_size)
...
queryPropertyTable(propertyTableLabel, jPropTable, rPropTable, cPropTable, jError)