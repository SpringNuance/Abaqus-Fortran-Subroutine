# hydra-utils

## Contents

`gencoh.py` is a python script for generating an Abaqus cohesive
material definition with a trapezoidal TSL.

## Installation and running

The code requires Python 2.7 or Python 3 and
[NumPy](http://www.numpy.org) 1.10 or newer.

        $ pip install --upgrade 'numpy>=1.10'

The scripts can be run without the need for installation, e. g.:

        $ python gencoh.py 0.00075 0.00975 0.0157 2600 > aba.inp
