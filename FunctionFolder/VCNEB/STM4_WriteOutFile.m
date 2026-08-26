function STM4_WriteOutFile(opt)
global POP_STRUC
global ORG_STRUC
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);
UnitTemp = 1.0;
if     0==opt
fpath = [ ORG_STRUC.resFolder '/PATH/path_' num2str(POP_STRUC.step) '.POSCAR'];
fp = fopen(fpath, 'w');
elseif 1==opt
fpath = [ ORG_STRUC.resFolder '/transitionPath_POSCARs' ];
fp = fopen(fpath, 'w');
elseif 2==opt
fpath = [ ORG_STRUC.resFolder, '/AuxiliaryFiles/gatheredPOSCARS'];;
fp = fopen(fpath, 'a');
end
for i = 1:numImages
if 2==opt
fprintf(fp, 'Step %3d -- Image  %2d \n', POP_STRUC.step, i);
else
fprintf(fp, 'Image  %2d \n', i);
end
fprintf(fp, '%12.10f\n', UnitTemp);
para = POP_STRUC.POPULATION(i).LATTICE/UnitTemp;
for j = 1:3
fprintf(fp, '    %16.12f   %16.12f  %16.12f\n', para(j,:));
end
fprintf(fp, ' %d ', ORG_STRUC.numIons(1:end) );
fprintf(fp, ' \nDirect \n');
for j = 1:sumIons
fprintf(fp, ' %16.12f   %16.12f  %16.12f\n', (POP_STRUC.POPULATION(i).COORDINATES(j,:)) );
end
end
fclose(fp);
