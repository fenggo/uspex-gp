function NEB_writeOUT_POSCAR(Ind_No)
global POP_STRUC
global ORG_STRUC
numIons = ORG_STRUC.numIons;
lat = POP_STRUC.POPULATION(Ind_No).LATTICE;
coord = POP_STRUC.POPULATION(Ind_No).COORDINATES;
bodyCount = POP_STRUC.bodyCount;
[nothing, nothing] = unix(['cat /dev/null > POSCAR' ]);
Lattice_par = latConverter(lat);
if size(Lattice_par,1) == 6
Lattice_par = Lattice_par';
end
L = sum(ORG_STRUC.numIons);
if size(coord,1) ~= L
error = 'number of variables (coordinates) of candidate structure do not correspond to the number of ions specified';
save ([ ORG_STRUC.homePath, '/', ORG_STRUC.resFolder '/ERROR_numVars.txt'],'error');
end
Lattice_par(4:6) = Lattice_par(4:6)*180/pi;
fpath = [ './POSCAR' ];
fp = fopen(fpath, 'w');
fprintf(fp, 'Image%d',Ind_No);
fprintf(fp, ' %6.4f ', Lattice_par(1:6));
fprintf(fp, '\n%3.10f\n', 1.0);
for latticeLoop = 1:3
fprintf(fp, ' %17.12f  %17.12f  %17.12f\n', lat(latticeLoop,:) );
end
fprintf(fp, ' %d ', numIons );
fprintf(fp, '\n');
fprintf(fp, 'Direct \n');
for coordLoop = 1:L
fprintf(fp, ' %17.12f  %17.12f  %17.12f\n', coord(coordLoop,:) );
end
fclose(fp);
return
