from __future__ import print_function
import sys
import os
import odbAccess as oa

def main():

    # parse command line arguments
    if not len(sys.argv) == 2:
        print('usage:\n    {0} <job name>'.format(sys.argv[0]),
              file=sys.stderr)
        sys.exit(2)
    odbpath = sys.argv[1] + '.odb'
    mappath = sys.argv[1] + '.map'

    print('%s -> %s' % (odbpath, mappath, ))

    # open odb
    try:
        odb = oa.openOdb(odbpath, readOnly=True)
    except oa.OdbError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

    # look for required element sets
    assembly = odb.rootAssembly
    try:
        instance = assembly.instances['PART-1-1']
    except KeyError as exc:
        print('odb does not have expected instance: {0}'.format(exc),
              file=sys.stderr)
        exit(1)
    try:
        map_c = instance.elementSets['MAP_C'].elements
        map_s = instance.elementSets['MAP_S'].elements
    except KeyError as exc:
        print('odb does not have expected element set: {0}'.format(exc),
              file=sys.stderr)
        exit(1)

    # write mapfile
    with open(mappath, 'w') as mapfile:
        print('## hydra generated map file ##', file=mapfile)
        for i, j in zip(map_s, map_c):
            print(i.label, j.label, file=mapfile)

if __name__ == '__main__':
    main()
