       PARAMETER(MWCOMP=number of wave components)
       DIMENSION WAMP(MWCOMP),WPERD(MWCOMP),WXLAMB(MWCOMP),
     1 WPHI(MWCOMP),WOFF(3),WANG(2,MWCOMP)
       ...
       CALL GETWAVE(MWCOMP,NWCOMP,WAMP,WPERD,WXLAMB,WPHI,WOFF,WANG,
     1 ELEVB,ELEVS,JWTYPE,JRCD)    


           CALL GETWAVEVEL (NDIM, X, V, A, LERROR, NOEL, XINTERMED)
      CALL GETWINDVEL (NDIM, X, V, NOEL, XINTERMED)
      CALL GETCURRVEL (NDIM, X, V, NOEL, XINTERMED)