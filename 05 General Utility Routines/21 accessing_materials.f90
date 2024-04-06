     dimension stress(*),ddsdde(ntens,*),stran(*),dstran(*),
    *	defGrad(3,3),predef(npredf),dpredef(npredf),coords(3)
     ...

     call material_lib_mech(materiallib,stress,ddsdde,stran,dstran,
    *       npt,dvdv0,dvmat,dfgrd,predef,dpredef,npredf,celent,coords)
     ...