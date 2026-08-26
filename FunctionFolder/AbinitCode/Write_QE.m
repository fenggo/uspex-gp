function Write_QE(Ind_No)
global POP_STRUC
global ORG_STRUC
numIons     = POP_STRUC.POPULATION(Ind_No).numIons;
Step        = POP_STRUC.POPULATION(Ind_No).Step;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;
try
[nothing, nothing] = unix(['cat qEspresso_options_' num2str(Step) ' > qe.in']);
[nothing, nothing] = unix(['sed -e "s/AAAA/' num2str(sum(numIons)) '/" qe.in > TEMP']);
[nothing, nothing] = unix(['sed -e "s/BBBB/' num2str(length(numIons)) '/" TEMP > qe.in']);
[nothing, nothing] = unix(['rm TEMP']);
catch
error = ['qEspresso_options is not present for step ' num2str(Step)];
quit
end
bohr = 0.529177; 
fp = fopen('qe.in','a+');
fprintf(fp, 'CELL_PARAMETERS cubic\n');
lat1 = latConverter(LATTICE);
lat = latConverter(lat1);
lat = lat/bohr;
fprintf(fp, '%8.4f %8.4f %8.4f\n', lat(1,:));
fprintf(fp, '%8.4f %8.4f %8.4f\n', lat(2,:));
fprintf(fp, '%8.4f %8.4f %8.4f\n', lat(3,:));
fprintf(fp, 'ATOMIC_POSITIONS {crystal} \n');
coordLoop = 1;
for i = 1 : length(numIons)
for j = 1 : numIons(i)
fprintf(fp, '%4s   %12.6f %12.6f %12.6f\n', megaDoof(ORG_STRUC.atomType(i)), COORDINATES(coordLoop,:));
coordLoop = coordLoop + 1;
end
end
[Kpoints, Error] = Kgrid(LATTICE, ORG_STRUC.Kresol(Step), ORG_STRUC.dimension);
if Error == 1  
POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
else
POP_STRUC.POPULATION(Ind_No).K_POINTS(Step,:)=Kpoints;
end
fprintf(fp, 'K_POINTS {automatic} \n');
fprintf(fp, '%4d %4d %4d  0 0 0\n', Kpoints(1,:));
fclose(fp);
