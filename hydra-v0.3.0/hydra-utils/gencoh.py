"""
generate material block cards for standard Hydra coehesive elements

example usage:

gencoh.py --help
gencoh.py 0.00075 0.00975 0.0157 2600 > COH.inp

--
Copyright (c) 2016 PoliHydra

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from __future__ import print_function
import numpy as np

COMMENT = {
    'START': """\
**\n** HYDRA GENERATED
** material   : {mat}
** delta_0    : {d0}
** delta_1    : {d1}
** delta_F    : {df}
** t0         : {t0}
** n01        : {n01}
** n1f        : {n1f}
** hydr. dec. : {hydrogen}
** tr. stress : {true_stress}
**""",
    }
KEYW = {
    'MATERIAL': '*MATERIAL, NAME={}',
    'ELASTIC': '*ELASTIC, TYPE=TRACTION, DEPENDENCIES={deps:d}',
    'DAMAGE INITIATION': '*DAMAGE INITIATION, CRITERION=MAXE',
    'DAMAGE EVOLUTION':
        '*DAMAGE EVOLUTION, TYPE=DISPLACEMENT, SOFTENING=TABULAR',
    'UFIELD': '*USER DEFINED FIELD',
    'DEPVAR': '*DEPVAR',
    }
DATA = {
    'ELASTIC': '{e_nn:G}, {e_ss:G}',
    'ELASTIC_DP':
        ('{e_nn:G}, {e_ss:G}',
         '{e_nn:.5E}, {e_ss:.5E},  ,  , {dep: .2f}',
         '{e_nn:.5E}, {e_ss:.5E},  ,  , {dep[0]: .2f}, {dep[1]: .2f}',
         ),
    'DAMAGE INITIATION': '{d0_n:G}, {d0_s:G}',
    'DAMAGE EVOLUTION': '{damage:f}, {delta:G}',
    }


def trapezoidal(d0, d1, df, n01=10, n1f=10):

    dam1 = 1 - d0/d1

    #
    # first compute dam from n01/2 linearly spaced delta
    #
    d01_l = np.linspace(d0, d1, n01//2, endpoint=False, )
    dam01 = 1 - d0/d01_l

    # check 0 <= dam01[i] < d1, with equality only for i==0
    assert dam01[0] == 0
    assert (0 < dam01[1:]).all() and (dam01 < dam1).all()

    #
    # then compute delta from n01/2 linearly spaced dam
    #
    dam01_l = np.linspace(0., dam1, n01//2, endpoint=False, )
    d01 = d0 / (1 - dam01_l)

    # check d0 <= d01[i] < d1, with equality only for i==0
    assert d01[0] == d0
    assert (d0 < d01[1:]).all() and (d01 < d1).all()

    # first point shold be the same for both computations
    assert (d01_l[0], dam01[0]) == (d01[0], dam01_l[0])

    #
    # finally merge points
    #
    d01 = np.concatenate((d01_l, d01[1:]))
    dam01 = np.concatenate((dam01, dam01_l[1:]))
    d01.sort()
    dam01.sort()

    # there can be duplicates/close duplicates
    dups = (np.isclose(d01[1:], d01[:-1]) | np.isclose(dam01[1:], dam01[:-1]))
    if dups.any():
        dupidx, = dups.nonzero()
        d01 = np.delete(d01, dupidx)
        dam01 = np.delete(dam01, dupidx)
    assert not np.isclose(d01[1:], d01[:-1]).any()
    assert not np.isclose(dam01[1:], dam01[:-1]).any()
    # sure to have done the math correctly?
    np.testing.assert_allclose(d01, d0 / (1 - dam01))
    np.testing.assert_allclose(dam01, 1 - d0/d01)

    #
    # descending curve
    #
    d1f = np.linspace(d1, df, n1f+1, )
    dam1f = 1 - d0/(df - d1)*(df/d1f - 1)
    assert (dam1f[:-1] < 1.0).all()
    assert dam1f[-1] == 1

    return np.hstack((d01, d1f, )), np.hstack((dam01, dam1f, ))


def geninp(d0, d1, df, t0, hydrogen=False, true_stress=False, mat='COH',
           n01=200, n1f=200, ):

    if not 0 < d0 < d1 < df:
        raise ValueError('d0, d1, df, values must be positive and increasing')
    if not 0 < t0:
        raise ValueError('t0 must be positive')

    ymodulus = t0/d0
    delta, dam = trapezoidal(d0, d1, df, n01=n01, n1f=n1f)

    print(COMMENT['START'].format(**vars()))
    print(KEYW['MATERIAL'].format(mat))

    dependencies = []
    if hydrogen:
        k = np.linspace(0.1, 1., num=2)
        factor = np.stack((k, k), axis=-1)
        dependencies.append(('K', factor))
    if true_stress:
        logE = np.linspace(-1., 1., num=51)
        factor = np.stack((logE, np.exp(logE)), axis=-1)
        dependencies.append(('LE', factor))

    if dependencies:
        print(KEYW['UFIELD'].format())
        print(KEYW['DEPVAR'].format())
        # fixme: use DATA dict
        print(len(dependencies))
        for i, (label, _) in enumerate(dependencies, 1):
            print(i, label, sep=', ')

    assert 0 <= len(dependencies) <= 2
    print(KEYW['ELASTIC'].format(deps=len(dependencies)))
    if len(dependencies) == 0:
        print(DATA['ELASTIC'].format(e_nn=ymodulus, e_ss=ymodulus))
    elif len(dependencies) >= 1:
        for ivar, ifactor in dependencies[-1][1]:
            yp = ymodulus*ifactor
            if len(dependencies) == 2:
                for jvar, jfactor in dependencies[-2][1]:
                    jyp = yp*jfactor
                    print(DATA['ELASTIC_DP'][2].format(
                        e_nn=jyp, e_ss=jyp, dep=[jvar, ivar]))
            else:
                print(DATA['ELASTIC_DP'][1].format(e_nn=yp, e_ss=yp, dep=ivar))

    print(KEYW['DAMAGE INITIATION'].format())
    print(DATA['DAMAGE INITIATION'].format(d0_n=d0, d0_s=d0))

    print(KEYW['DAMAGE EVOLUTION'].format())
    for i, j in zip(dam, delta-d0):
        print(DATA['DAMAGE EVOLUTION'].format(damage=i, delta=j))


def main():
    """main entry point"""
    import argparse

    parser = argparse.ArgumentParser(description='generate COH material block')
    parser.add_argument('d0', metavar='delta_0', type=float, )
    parser.add_argument('d1', metavar='delta_1', type=float, )
    parser.add_argument('df', metavar='delta_F', type=float, )
    parser.add_argument('t0', metavar='sigma_0', type=float, )
    parser.add_argument('--hydrogen', action='store_true',
                        help='apply K factor dependency'
                             ' (default no K dependency)')
    parser.add_argument('--true-stress', action='store_true',
                        help='apply true stress correction to TSL'
                             ' (default use nominal stress)')

    args = parser.parse_args()
    geninp(**vars(args))


if __name__ == '__main__':
    main()
