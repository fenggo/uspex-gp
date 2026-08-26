import json
import os
import spglib
import sys

ERROR = 0
N_atoms = 0
if __name__ == '__main__':
    structurePath = str(sys.argv[1])

    assert os.path.exists(structurePath)
    with open(structurePath, 'r') as fp:
        tmp = json.loads(fp.readline())

    if isinstance(tmp['numbers'], int):
        tmp['numbers'] = [tmp['numbers']]
        tmp['coordinates'] = [tmp['coordinates']]
    cell = (tmp['lattice'], tmp['coordinates'], tmp['numbers'])
    try:
        new_cell, new_frac_positions, new_numbers = spglib.standardize_cell(cell, to_primitive=True, symprec=tmp['tolerance'])
    except TypeError:
        new_cell, new_frac_positions, new_numbers = spglib.standardize_cell(cell, to_primitive=True)
    N_atoms = len(new_numbers)
    try:
        from ase import Atoms
        from ase.io.cif import write_cif
        from ase.io.vasp import write_vasp
        new_atoms = Atoms(numbers=new_numbers, cell=new_cell, scaled_positions=new_frac_positions)
        write_cif('symmetrized_spg.cif', new_atoms)
        write_vasp('outputPOSCAR', new_atoms, direct=True, vasp5=True)

    except ImportError:
        ERROR = 1
        print('Cannot import ase. Spglib symmetrization will be skipped! Stokes symmetrization will be used!')

    print('<CALLRESULT> ' + str(N_atoms) + ' ' + str(ERROR))
