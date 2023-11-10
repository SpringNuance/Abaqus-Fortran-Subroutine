from __future__ import print_function
import sys
import os
import odbAccess as oa
import numpy as np

HYDRA = 'HYDRA'

def main():

    if not len(sys.argv) == 2:
        print('usage:\n    {0} <job name>'.format(sys.argv[0]),
              file=sys.stderr)
        sys.exit(2)
    odbpath = sys.argv[1] + '.odb'
    jacpath = sys.argv[1] + '.jac'

    print('%s -> %s' % (odbpath, jacpath,))

    odb = oa.openOdb(odbpath, readOnly=True)
    instances = odb.rootAssembly.instances

    b = dict()
    b['CPx4R'] = 0.25 * np.array(((-1., +1., +1., -1.),
                                  (-1., -1., +1., +1.),),)

    with open(jacpath, 'w') as txt:
        txt.write('** inverse jacobian matrix\n')
        for element in instances['PART-1-1'].elementSets[HYDRA].elements:
                if element.type in ['CPS4R', 'CPE4R', 'CPE4RT']:
                    jac = np.zeros((2,2))
                    xy = np.zeros((4,2))
                    inames = element.instanceNames
                    conn = element.connectivity
                    for n in range(4):
                        xy[n,:] = instances[inames[n]].getNodeFromLabel(
                                conn[n]).coordinates[0:2]
                    jac[:,:] -= xy[0,:]
                    jac[0,:] += xy[1,:]
                    jac[1,:] -= xy[1,:]
                    jac[:,:] += xy[2,:]
                    jac[0,:] -= xy[3,:]
                    jac[1,:] += xy[3,:]
                    jac *= 0.25
                    grad = np.linalg.solve(jac, b['CPx4R'])
                    txt.write("** {0}\n".format(element.type))
                    txt.write("%d\n" % element.label)
                    np.savetxt(txt, conn, fmt='%d')
                    np.savetxt(txt, grad)
                else:
                    raise NotImplementedError, element.type

if __name__ == '__main__':
    main()
