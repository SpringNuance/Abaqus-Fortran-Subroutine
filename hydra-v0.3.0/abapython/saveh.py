from __future__ import print_function
import sys

import odbAccess as oa
import numpy as np
import numpy.testing
import re

elip = re.compile(r'Element PART-1-1.(?P<el>\d+) Int Point (?P<ip>\d+)').match
node = re.compile(r'Node PART-1-1.(?P<nd>\d+)').match

labelmap = [(elip, lambda x: 'E{el}IP{ip}_'.format(**x.groupdict())),
            (node, lambda x: 'N{nd}_'.format(**x.groupdict())), ]


def genlab(reg):
    for matf, genl in labelmap:
        match = matf(reg)
        if match:
            return genl(match)
    else:
            return None


def main():

    path = sys.argv[1]

    dbpath = path + '.odb'
    savepath = path + '.npz'

    odb = oa.openOdb(dbpath, readOnly=True)

    print('-'*60)
    print('name:', odb.name)
    print('analysisTitle:', odb.analysisTitle)
    print('creationTime:', odb.jobData.creationTime)

    step = odb.steps['Step-1']

    print('step name:', step.name)
    print('step procedure:', step.procedure)
    print('step description:', step.description)
    print('step domain:', step.domain)
    print('step frames: %d' % (len(step.frames), ))

    t = np.fromiter((i.frameValue for i in step.frames), dtype=np.float)

    res = {'t': t}
    for hr in step.historyRegions.keys():
        label = genlab(hr)
        if not label:
            continue
        houts = step.historyRegions[hr].historyOutputs
        for key in houts.keys():
            t1, v = np.asarray(houts[key].data).T
            np.testing.assert_array_equal(t, t1)
            res[label+key] = v

    print('saving %s to %s' % (sorted(res.keys()), savepath))
    np.savez(savepath, **res)
    print('-'*60)

if __name__ == '__main__':
    main()
