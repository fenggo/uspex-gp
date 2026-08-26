function META_extractStructures(max_num, weight)
global USPEX_STRUC
atomType = USPEX_STRUC.atomType;
N = length(USPEX_STRUC.POPULATION);
count = 1;
for i=1:N
if USPEX_STRUC.POPULATION(i).ToCount > 0
ID(count)   = i;
numIons     = USPEX_STRUC.POPULATION(i).numIons;
enth(count) = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(numIons);
count = count + 1;
end
end
if ~exist('enth')
disp('ERROR! Please check your INCAR and/or POSCAR files. Exit.');
quit();
end
[nothing, ranking] = sort(enth);
if count > 2
for i = 2 : count-1
if USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount > 0
for j= 1 : i-1
if USPEX_STRUC.POPULATION(ID(ranking(j))).ToCount > 0
if SameStructure(ID(ranking(i)), ID(ranking(j)), USPEX_STRUC)
USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount = 0;
break;
end
end
end
end
end
end
GoodList = zeros(max_num, 1); 
item = 1;
GoodList(1) = 1; 
if count > 2
for i = 2 : count-1
if USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount > 0
same = 0;
for j= 1 : i-1
if USPEX_STRUC.POPULATION(ID(ranking(j))).ToCount > 0
if SameStructure_order(ID(ranking(i)), ID(ranking(j)), USPEX_STRUC)
USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount = 0;
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
end
Good_num = item; 
fp1 = fopen('goodStructures', 'w');
fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes     SYMM\n');
fprintf(fp1,'                       (eV/atom)    (A^3/atom)       \n');
[nothing, nothing] = unix(['cat /dev/null > goodStructures_POSCARS']);  
for i = 1 : Good_num
j =  ID(ranking(GoodList(i)));
num    = USPEX_STRUC.POPULATION(j).numIons;
volume = det(USPEX_STRUC.POPULATION(j).LATTICE)/sum(num);
symg   = USPEX_STRUC.POPULATION(j).symg;
lattice= USPEX_STRUC.POPULATION(j).LATTICE;
coor   = USPEX_STRUC.POPULATION(j).COORDINATES;
composition = sprintf('%3d',num);
shift=[4, 2, 1]; 
if size(composition,2)<11
composition=[composition,blanks(shift(length(num)))];
end
fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %4d\n',...
j, composition, enth(ranking(GoodList(i))), volume, symg);
Write_POSCAR(atomType, j, symg, num, lattice, coor);
[nothing, nothing] = unix([' cat POSCAR      >> goodStructures_POSCARS']);
end
fclose(fp1);
