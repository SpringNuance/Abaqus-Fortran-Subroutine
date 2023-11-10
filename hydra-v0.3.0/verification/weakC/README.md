# Weak coupling verification problem

## 3-step procedure

### Abaqus job files:

* H-1static.inp: preliminary static analysis job file
* H-2massdiff.inp: mass diffusion analysis job file
* H-3cohesive.inp: cohesive analysis job file
* _AISI.inp: include file with AISI steel material block
* _cohes.inp: include file with cohesive mesh
* _equation.inp: include file with *EQUATION cards
* _geometry.inp: include file with node geometry
* _solid.inp: continuum elements mesh for steps 1static and 3cohesive
* _solid-2massdiff.inp: continuum elements mesh for mass diffusion analysis
* conc.f: user subroutine
* conc.inc: user subroutine include file

### Hydra generated files:

* _coh-H.inp: include file with cohesive TSL law (generated)
* H-3cohesive.map: solid to cohesive mapping file (generated)

### Run instructions.

1. Perform a preliminary static analysis to calculate the hydrostatic
   stress field.
   In this example a gravity load is applied in order to obtain a vertical
   stress gradient over the elements

        $ abaqus analysis job=H-1static

2. Mass diffusion analysis

        $ abaqus analysis job=H-2massdiff

3. Generate cohesive element TSL

        $ python ../../hydra-utils/gencoh.py 0.22e-3 9.7e-3 15.e-3 1000 --hydrogen > _coh-H.inp

4. Generate the solid-to-cohesive map file (`H-3cohesive.map`)

        $ abaqus datacheck job=H-3cohesive
        $ abaqus python ../../abapython/mkmap.py H-3cohesive

5. Perform the cohesive analysis

        $ abaqus continue job=H-3cohesive user=conc.f

6. Post-process results

        $ abaqus python ../../abapython/saveh.py H-3cohesive

Reference values for the expected results can be found in the `results`
directory. Initial normalized concentration (`NNC`) can be changed in the
preamble of file `H-2massdiff.inp`, under the `*PARAMETER` card.

Results for the following initial concentrations are reported:

* `NNC=7.042` in `results/H-3cohesive_NNC=07.npz`
* `NNC=14.084` in `results/H-3cohesive_NNC=14.npz`
* `NNC=21.127` in `results/H-3cohesive_NNC=21.npz`
