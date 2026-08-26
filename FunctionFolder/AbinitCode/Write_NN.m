function writeCalcFilesNN(Ind_No)
global POP_STRUC
global ORG_STRUC
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
L = sum(numIons);
if size(POP_STRUC.POPULATION(Ind_No).COORDINATES,1) ~= L
error = 'number of variables (coordinates) of candidate structure do not correspond to the number of ions specified';
save ([ ORG_STRUC.resFolder '/ERROR_numVars.txt'],'error')
end
bohr = 0.529177; 
lat = POP_STRUC.POPULATION(Ind_No).LATTICE/bohr;  
[nothing, nothing] = unix(['echo ' num2str(lat(1,:)) ' > structure.dat']);
[nothing, nothing] = unix(['echo ' num2str(lat(2,:)) ' >> structure.dat']);
[nothing, nothing] = unix(['echo ' num2str(lat(3,:)) ' >> structure.dat']);
for i = 1:L
[nothing, nothing] = unix(['echo ' num2str(POP_STRUC.POPULATION(Ind_No).COORDINATES(i,:)*lat) ' >> structure.dat']);
end
