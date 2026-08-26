function writeOUT_POSCAR(Ind_No, doOrder)
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
Lattice_par(4:6) = Lattice_par(4:6)*180/pi;
[nothing, nothing] = unix(['echo ' 'EA' num2str(bodyCount) ' ' num2str(Lattice_par) '  >> POSCAR']);
[nothing, nothing] = unix(['echo ' '1.0' ' >> POSCAR']);
L = sum(ORG_STRUC.numIons);
if size(coord,1) ~= L
error = 'number of variables (coordinates) of candidate structure do not correspond to the number of ions specified';
save ([ ORG_STRUC.resFolder '/ERROR_numVars.txt'],'error');
end
for latticeLoop = 1 : 3
[nothing, nothing] = unix(['echo ' num2str(lat(latticeLoop,:),12) ' >> POSCAR']);
end
[nothing, nothing] = unix(['echo ' num2str(numIons) ' >> POSCAR']);
[nothing, nothing] = unix(['echo ' 'Direct' ' >> POSCAR']);
for coordLoop = 1 : L
[nothing, nothing] = unix(['echo ' num2str(coord(coordLoop,:),12) ' >> POSCAR']);
end
