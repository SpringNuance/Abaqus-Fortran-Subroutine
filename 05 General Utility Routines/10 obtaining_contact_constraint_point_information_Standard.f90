PARAMETER(NDATASIZE=10)
DIMENSION ARRAY(NDATASIZE), JARRAY(NDATASIZE)
CHARACTER*3 FLGRAY(NDATASIZE)
...
CALL GETVRC('VAR', ARRAY, JARRAY, FLGRAY, NDATASIZE, JRCD)