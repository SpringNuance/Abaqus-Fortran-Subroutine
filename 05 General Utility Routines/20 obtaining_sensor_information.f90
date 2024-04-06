     character*80 mySensorName
     ...

     iMySensorID =  IGETSENSORID(mySensorName, jSensorLookUpTable )
     iMySensorID = IVGETSENSORID(mySensorName, jSensorLookUpTable )
     dMySensorValue = sensorValues(iMySensorID)
     ...
     dMySensorValue =   GETSENSORVALUE(mySensorName, 
C                       jSensorLookUpTable, sensorValues )
     dMySensorValue =  VGETSENSORVALUE(mySensorName,  
C                      jSensorLookUpTable, sensorValues )
     ...