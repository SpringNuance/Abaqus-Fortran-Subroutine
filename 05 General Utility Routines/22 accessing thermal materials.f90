     dimension predef(npredef),dpredef(npredef),dtemdx(*),
    *         rhodUdg(*),flux(*),dfdt(*),dfdg(ndim,*),drpldt(*),    *         coords(3)
     ...

     call material_lib_ht(materiallib,rhoUdot,rhodUdt,rhodUdg,
    *       flux,dfdt,dfdg,rpl,drpldt,npt,dvmat,predef,
    *       dpredef,npredf,temp,dtemp,dtemdx,celent,coords)
     ...