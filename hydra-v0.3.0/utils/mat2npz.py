from __future__ import print_function
import os
import argparse
from numpy import savez as savenpz
from scipy.io import loadmat


def main():
    parser = argparse.ArgumentParser(
        description='convert matlab .mat to numpy .npz format')
    parser.add_argument('mat', nargs='+', metavar='input.mat',
                        help='input .mat file')

    args = parser.parse_args()

    for mat in args.mat:
        root, _ = os.path.splitext(mat)
        npz = root + '.npz'
        try:
            data = loadmat(mat, squeeze_me=True,)
            savenpz(npz, **data, )
        except (IOError, NotImplementedError) as exp:
            print('%s: %s' % (mat, exp))
        else:
            print('%s -> %s' % (mat, npz))


if __name__ == '__main__':
    main()
