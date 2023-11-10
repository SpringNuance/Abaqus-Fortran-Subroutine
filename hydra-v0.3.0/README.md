# hydra

**HYDR**ogen embrittlement in **A**baqus FEA: a collection of tools for
the simulation of hydrogen embrittlement by means cohesive elements and
Abaqus user subroutines.


## Directory structure

* `abafortran`: *hydra* Abaqus user subroutines
* `abapython`: *hydra* Abaqus python utilities
* `hydra-utils`: *hydra* generic python utilities
* `verification`: verification examples for testing *hydra*
* `examples`: example problems illustrating *hydra* capabilities
* `utils`: useful python utilities, not directly linked to *hydra*

# Installing and running the examples

*hydra* is a collection of utilities which can be used without the
need of a proper installation, however some prerequisites are needed:

- an
  [Abaqus/Standard](https://www.3ds.com/products-services/simulia/products/abaqus/abaqusstandard/)
  installation with a supported FORTRAN compiler. The present *hydra*
  version has been tested with Abaqus versions from 6.12-1 to 2017.
- a fairly recent [python](https://www.python.org) installation, in addition
  to the python version shipped with Abaqus. Currently version 2.7 or
  3.6 are suggested.
  
For the example problems in the `examples` and `verification`
directories a `README` file is present, with step-by-step instructions
for executing the described analyses. These instructions are written
under the assumption that all examples are run within the example
directory itself.

## Compatibility issues

- symlinks: as *hydra* was developed under a Unix-like (actually
  mostly Posix compliant) operating system, its source distribution
  contains [symbolic
  links](https://en.wikipedia.org/wiki/Symbolic_link). Symlinks are
  fully supported only in Windows 10, see
  <https://blogs.windows.com/buildingapps/2016/12/02/symlinks-windows-10/>.
  Under previous versions of Windows, symbolic links may be get
  substituted by a small text file containing the name of the file
  they are pointing to. In such cases simply copy the original file in
  place of the symlink.
