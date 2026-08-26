function extractStructures(max_num, weight)
global USPEX_STRUC
atomType = USPEX_STRUC.SYSTEM(1).atomType;
N = length(USPEX_STRUC.POPULATION);
for i=1:N
numIons    = USPEX_STRUC.POPULATION(i).numIons;
fitness(i) = USPEX_STRUC.POPULATION(i).Fitness;
enth(i)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(numIons);
end
[nothing, ranking] = sort(fitness);
GoodList = zeros(max_num, 1); 
item = 1;
GoodList(1) = 1; 
for i = 2 : N
if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
same = 0;
for j= 1 : i-1
if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
if SameStructure(ranking(i), ranking(j),USPEX_STRUC)
USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
same = 1;
break;
end
end
end
if same == 0
item = item + 1;
GoodList(item) = i;
end
end
if item == max_num
break;
end
end
Good_num = item; 
fp1 = fopen('goodStructures', 'w');
[a,b]=unix(['cat /dev/null > goodStructures_POSCARS']);  
fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes    fitness  SYMM\n');
fprintf(fp1,'                       (eV/atom)    (A^3/atom)    ()         \n');
for i = 1 : Good_num
j = ranking(GoodList(i));
num    = USPEX_STRUC.POPULATION(j).numIons;
volume = USPEX_STRUC.POPULATION(j).Vol/sum(num);
symg   = USPEX_STRUC.POPULATION(j).symg;
lattice= USPEX_STRUC.POPULATION(j).LATTICE;
coor   = USPEX_STRUC.POPULATION(j).COORDINATES;
composition = sprintf('%3d',num);
shift=[4, 2, 1]; 
if size(composition,2)<11
composition=[composition,blanks(shift(length(num)))];
end
fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %9.4f  %4d\n',...
j, composition, enth(j), volume, fitness(j), symg);
Write_POSCAR(atomType, j, symg, num, lattice, coor);
[nothing, nothing] = unix([' cat POSCAR      >> goodStructures_POSCARS']);
end
fclose(fp1);
