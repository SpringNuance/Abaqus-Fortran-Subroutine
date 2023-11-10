from __future__ import print_function
import os
import argparse
from numpy import load
from scipy.io import savemat


def main():
    parser = argparse.ArgumentParser(
        description='convert numpy .npz to matlab .mat format')
    parser.add_argument('npz', nargs='+', metavar='input.npz',
                        help='input .npz file')

    args = parser.parse_args()

    for npz in args.npz:
        root, _ = os.path.splitext(npz)
        mat = root + '.mat'
        try:
            with load(npz) as data:
                savemat(mat, data, )
        except IOError as exp:
            print(exp)
        else:
            print('%s -> %s' % (npz, mat))


if __name__ == '__main__':
    main()
