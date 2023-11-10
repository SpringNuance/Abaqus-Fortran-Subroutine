import sys
import argparse
import os.path
import numpy as np
import numpy.linalg


def norm(x):
    return numpy.linalg.norm(x, ord=np.inf)


def main():
    parser = argparse.ArgumentParser(
        description='compare a numpy .npz file to a reference file')
    parser.add_argument('data', metavar='A.npz', help='file to test')
    parser.add_argument('refv', metavar='B.npz', help='reference file')

    args = parser.parse_args()

    if os.path.isdir(args.refv):
        refv = os.path.join(args.refv, os.path.basename(args.data))
    else:
        refv = args.refv
    try:
        with np.load(args.data) as data, np.load(refv) as ref:
            dataset = set(data.files)
            refset = set(ref.files)
            equalset = set()
            for key in sorted(refset & dataset):
                a = data[key]
                b = ref[key]
                if a.shape != b.shape:
                    print('* {}: different shape {} {}.'.format(
                        key, a.shape, b.shape))
                    continue
                if np.all(a == b):
                    equalset.add(key)
                else:
                    try:
                        err = norm(a-b)
                    except TypeError:
                        print('* {}: differ'.format(key))
                    else:
                        print('* {}: |diff|_inf = {}'.format(key, err))
            if refset - dataset:
                print('* missing vars in {}:'.format(args.data))
                for k in sorted(refset - dataset):
                    print('   {}'.format(k))
            if dataset - refset:
                print('* extra vars in {}:'.format(args.data))
                for k in sorted(dataset - refset):
                    print('   {}'.format(k))
            if equalset == dataset:
                print("Files '{}' and '{}' are identical.".format(
                    args.data, refv))
                sys.exit(0)
            else:
                print('* other {} vars identical.'.format(len(equalset)))
                sys.exit(1)
    except IOError as exp:
        print(exp)


if __name__ == '__main__':
    main()
