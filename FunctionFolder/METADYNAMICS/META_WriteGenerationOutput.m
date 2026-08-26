function META_WriteGenerationOutput(fitness, flag)
global POP_STRUC
global ORG_STRUC
global USPEX_STRUC
atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;
IND=POP_STRUC.ranking(1);
lattice = POP_STRUC.POPULATION(IND).LATTICE;
coor    = POP_STRUC.POPULATION(IND).COORDINATES;
numIons = POP_STRUC.POPULATION(IND).numIons;
lattice = POP_STRUC.POPULATION(IND).LATTICE;
symg    = POP_STRUC.POPULATION(IND).symg;
order   = POP_STRUC.POPULATION(IND).order;
count   = POP_STRUC.POPULATION(IND).Number;
if flag == 1
Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
[nothing, nothing] = unix([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
[nothing, nothing] = unix([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);
else
ID = update_USPEX_GENERATION(IND, ORG_STRUC.varcomp, ORG_STRUC.conv_till);
META_WriteOUTPUT(ID, resFolder, ORG_STRUC.varcomp, ORG_STRUC.conv_till, 3);
if (ORG_STRUC.FullRelax > 0)
Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
[nothing, nothing] = unix([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS_relaxed']);
META_WriteOUTPUT(ID, resFolder, ORG_STRUC.varcomp, length([ORG_STRUC.abinitioCode]), 4);
end
safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
cd (resFolder);
META_makeFigures(ORG_STRUC.conv_till);
if (ORG_STRUC.FullRelax > 0)
META_extractStructures(ORG_STRUC.populationSize, ORG_STRUC.weight);
end
cd ..
end
