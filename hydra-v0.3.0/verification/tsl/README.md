# TSL verification problems

## noH

Single element test for verifying the correct implementation of a trapezoidal
TSL.

### Abaqus job files:

* `noH.inp`: job file
* `_AISI.inp`: include file with AISI steel material block
* `_mesh.inp`: include file with mesh definition

### Hydra generated files:

* `_coh.inp`: include file with cohesive TSL law (generated)

### Analysis result files:

* `noH.odb`: abaqus output database
* `noH.npz`: selected results in numpy format

### Run instructions.

1. Generate the `_coh.inp` material block using the `hydra-utils/gencoh.py` program

        $ python ../../hydra-utils/gencoh.py 0.22e-3 9.7e-3 15e-3 1000 > _coh.inp
   
   The  generated file should be identical to `./results/_coh.inp`

2. Run abaqus

        $ abaqus analysis job=noH

3. Post-process results

        $ abaqus python ../../abapython/saveh.py noH

   Post-processing will generate the data file `noH.npz`, which is stored
   in numpy [`.npz` format](http://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html).
   A reference results file is given in the `results` directory.


## noH-ts

Single element test for verifying the correct implementation of a trapezoidal
TSL with true stress (Cauchy) correction.

### Abaqus job files:

* `noH-ts.inp`: job file
* `_AISI.inp`: include file with AISI steel material block
* `_mesh.inp`: include file with mesh definition
* `ts.f`: user subroutine
* `ts.inc`: user subroutine include file

### Hydra generated files:

* `_coh-ts.inp`: include file with cohesive TSL law (generated)
* `noH-ts.map`: solid to cohesive mapping file (generated)

### Analysis result files:

* `noH-ts.odb`: abaqus output database
* `noH-ts.npz`: selected results in numpy format

### Run instructions.

1. Generate the `_coh-ts.inp` material block using the `hydra-utils/gencoh.py` program

        $ python ../../hydra-utils/gencoh.py 0.22e-3 9.7e-3 15e-3 1000 --true-stress > _coh-ts.inp

   The  generated file should be identical to `./results/_coh.inp`

2. Generate mesh `noH-ts.map` file

        $ abaqus datacheck j=noH-ts
        $ abaqus python ../../abapython/mkmap.py noH-ts

   The  generated file should be identical to `./results/noH-ts.map`

3. Run abaqus

        $ abaqus continue j=noH-ts user=ts.f

4. Post-process results

        $ abaqus python ../../abapython/saveh.py noH-ts

   Post-processing will generate the data file `noH-ts.npz`.
   A reference results file is given in the `results` directory.
