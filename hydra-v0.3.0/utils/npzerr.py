import argparse
import os.path
import numpy as np
import scipy as sp
import scipy.interpolate


def norm(x, t):
    dif = sp.integrate.cumtrapz(np.abs(x), t, initial=0.)[-1]
    return dif / (t[-1] - t[0])


def main():
    parser = argparse.ArgumentParser(
        description='compute relative error between a .npz file and a '
                    'reference file')
    parser.add_argument('data', metavar='A.npz', help='file to test')
    parser.add_argument('refv', metavar='B.npz', help='reference file')

    args = parser.parse_args()

    if os.path.isdir(args.refv):
        refvpath = os.path.join(args.refv, os.path.basename(args.data))
    else:
        refvpath = args.refv
    try:
        with np.load(args.data) as data, np.load(refvpath) as refv:
            dataset = set(data.files) - {'t'}
            refset = set(refv.files) - {'t'}
            if refset - dataset:
                print('* missing vars in {}:'.format(args.data))
                for k in sorted(refset - dataset):
                    print('   {}'.format(k))
            if dataset - refset:
                print('* extra vars in {}:'.format(args.data))
                for k in sorted(dataset - refset):
                    print('   {}'.format(k))
            # compute common time axis
            tdata = data['t']
            trefv = refv['t']
            t = np.union1d(tdata, trefv)
            # compute error
            equalset = set()
            for key in sorted(refset & dataset):
                # interpolate data on common time axis
                a = data[key]
                b = refv[key]
                aintp = sp.interpolate.interp1d(tdata, a, kind='linear')(t)
                bintp = sp.interpolate.interp1d(trefv, b, kind='linear')(t)
                # compute error norm
                err = norm(aintp-bintp, t)
                nor = norm(bintp, t)
                if err:
                    print('{:15} {:<8.6f} {:.1e}'.format(
                        key, err/nor, nor))
                else:
                    print('{:15} {:<8.0f} {:.1e}'.format(
                        key, err, nor))
                    equalset.add(key)
            if equalset == dataset:
                print("* files '{}' and '{}' are identical.".format(
                    args.data, refvpath))
    except IOError as exp:
        print(exp)

if __name__ == '__main__':
    main()
