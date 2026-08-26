function META_WriteIndividualOutput(Ind_No, Step, flag)
global POP_STRUC
global ORG_STRUC
varcomp   = ORG_STRUC.varcomp;
atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;
symg        =POP_STRUC.POPULATION(Ind_No).symg;
count       =POP_STRUC.POPULATION(Ind_No).Number;
lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
coor    = POP_STRUC.POPULATION(Ind_No).COORDINATES;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
order   = POP_STRUC.POPULATION(Ind_No).order;
if ORG_STRUC.doSpaceGroup == 1
[nothing, current_path] = unix('pwd');
current_path(end) = [];
cd([ORG_STRUC.homePath '/CalcFoldTemp']);
POP_STRUC.POPULATION(Ind_No).symg = anasym_stokes(lattice, coor, numIons, atomType, ORG_STRUC.SGtolerance);
symg        =POP_STRUC.POPULATION(Ind_No).symg;
cd(current_path)
[nothing, nothing] = unix(['echo data_findsym-STRUC-' num2str(count) '            >> ' resFolder '/symmetrized_structures.cif']);
[nothing, nothing] = unix(['cat CalcFoldTemp/symmetrized.cif                    >> ' resFolder '/symmetrized_structures.cif']);
end
Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
if flag == 1
[nothing, nothing] = unix([' cat POSCAR      >> ' POP_STRUC.resFolder '/gatheredPOSCARS']);
Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
[nothing, nothing] = unix([' cat POSCAR_order >>' POP_STRUC.resFolder '/gatheredPOSCARS_order']);
update_USPEX_INDIVIDUAL(POP_STRUC.POPULATION(Ind_No), resFolder, ...
POP_STRUC.generation, atomType, 1);
else 
[nothing, nothing] = unix([' cat POSCAR  >>     ' resFolder '/gatheredPOSCARS_relaxed']);
update_USPEX_INDIVIDUAL(POP_STRUC.POPULATION(Ind_No), resFolder, ...
POP_STRUC.generation, atomType, 2);
end
META_WriteOUTPUT(count, resFolder, varcomp, Step, flag);
