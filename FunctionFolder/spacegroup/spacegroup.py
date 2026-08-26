import json
import os
import spglib
import sys


if __name__ == '__main__':
    structurePath = str(sys.argv[1])

    assert os.path.exists(structurePath)
    with open(structurePath, 'r') as fp:
        tmp = json.loads(fp.readline())

    if isinstance(tmp['numbers'], int):
        tmp['numbers'] = [tmp['numbers']]
        tmp['coordinates'] = [tmp['coordinates']]
    cell = (tmp['lattice'], tmp['coordinates'], tmp['numbers'])
    spgName, spgNumber = spglib.get_spacegroup(cell, symprec=tmp['tolerance']).split()
    spgNumber = spgNumber.replace('(', '').replace(')', '')
    print('<CALLRESULT>')
    print(int(spgNumber))
