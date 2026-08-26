function extractStableStructures(ids)
global ORG_STRUC
global USPEX_STRUC
[nothing, nothing] = unix('rm -f POSCAR stablePOSCARS_pressure');
for i=1:length(ids)
Write_POSCAR(ORG_STRUC.atomType, ids(i),             ...
USPEX_STRUC.POPULATION(ids(i)).symg,    ...
USPEX_STRUC.POPULATION(ids(i)).numIons, ...
USPEX_STRUC.POPULATION(ids(i)).LATTICE, ...
USPEX_STRUC.POPULATION(ids(i)).COORDINATES);
[nothing, nothing] = unix('cat POSCAR >> stablePOSCARS_pressure');
end
[nothing, nothing] = unix('rm -f POSCAR');
end
