# C(T) specimen, fully coupled analysis.

In this example a fully-coupled analysis of a C(T) specimen is
presented.  This example is discussed in detail in

> **G. Gobbi, C. Colombo, S. Miccoli, L. Vergani (2018)**. A
> fully coupled implementation of hydrogen embrittlement in FE
> analysis. To be submitted.

## Files

### Abaqus Job Files

- `H.inp` analysis setup
- `_AISI.inp` material definition
- `_mesh.inp` mesh definition

### Abaqus user subrotines

In this example two different laws were adopted for the estimation of
the number of trap sites, Kumnick and Johnson (1980), here codenamed
`_k`, and Sofronis et. al. (2001), codenamed `_s`.

- `hydra_k.f`, UMATHT plus auxiliary routines for Kumnick law,
- `hydra_s.f`, same for Sofronis law,
- `hydra.inc`, include file.

### Hydra generated files

- `_coh.inp` TSL for cohesive elements
- `H.jac` inverse Jacobian matrix
- `H.map` continuum to cohesive elements mapping

## Run instructions

### Generation of mapping and Jacobian files

These steps are common for both trap sites laws.

1. Generate the TSL law for a coupled analysis:

        $ python ../../hydra-utils/gencoh.py 0.00075 0.00975 0.0157 2600 --hydrogen _coh.inp

   This command will generate `_coh.inp`

1. Generate the inverse Jacobian matrix and continuum to cohesive mapping

        $ abaqus j=H datacheck
		$ abaqus python ../../abapython/mkmap.py H
		$ abaqus python ../../abapython/mkjac.py H

   These commands will generate `H.jac` and `H.map`.

For reference the generated files are included in the
[results](results) directory.

### Coupled analysis

1. Run the coupled analysis with Kumnick law:

        $ ln H.jac H_k.jac
		$ ln H.map H_k.map
		$ abaqus j=H_k input=H user=hydra_k.f

1. same with Sofronis law:

		$ ln H.jac H_s.jac
		$ ln H.map H_s.map
		$ abaqus j=H_s input=H user=hydra_s.f



If a `***ERROR: MAXIMUM SIZE OF STATIC WORKSPACE HAS BEEN EXCEEDED.`
error occurs, please see the discussion for the weakly coupled
analysis [../CT](../CT)
