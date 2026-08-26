function extendedConvexHull_301(convex_hull, fit_max)
global POOL_STRUC
global ORG_STRUC
global USPEX_STRUC
atomType = ORG_STRUC.atomType;
stopCrit = ORG_STRUC.stopCrit;
popSize  = max(ORG_STRUC.populationSize, 20);  
bestHM   = max(ORG_STRUC.keepBestHM, 4);       
numIons  = ORG_STRUC.numIons;
weight   = ORG_STRUC.weight;
N        = length(USPEX_STRUC.POPULATION);
generation  = length(USPEX_STRUC.GENERATION);
[Nb, Ntype] = size(numIons);
enth    = zeros(N,1);
fitness = zeros(N,1);
comp    = zeros(N,Ntype);
comp1   = zeros(N,Nb);
for i=1:N
comp(i,:)  = USPEX_STRUC.POPULATION(i).numIons;    
comp1(i,:) = USPEX_STRUC.POPULATION(i).numBlocks;  
enth(i)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(comp(i,:));
end
fitness = final_convex_hull(convex_hull, enth, comp1, numIons);
for i=1:N
USPEX_STRUC.POPULATION(i).Fitness = fitness(i);
end
[nothing, ranking] = sort(fitness);
accepted = zeros(1, N);
newInd   = zeros(1, N);
Na = 0; 
for i = 1 : N
if fitness(ranking(i)) <= fit_max
Na = Na + 1;
if fitness(ranking(i)) >= -0.05 
accepted(ranking(i)) = 1;
if USPEX_STRUC.POPULATION(ranking(i)).gen == generation
newInd(ranking(i)) = 1;
end
end
else
break;
end
end
for i = 2 : Na
if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
for j = 1 : i-1
if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
if sameComposition(comp(ranking(i),:),comp(ranking(j),:))
if SameStructure(ranking(i), ranking(j), USPEX_STRUC)
USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
break;
end
end
end
end
end
end
for i = 2 : Na
if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
for j = 1 : i-1
if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
if sameComposition(comp(ranking(i),:),comp(ranking(j),:))
if SameStructure_order(ranking(i), ranking(j), USPEX_STRUC)
USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
accepted(ranking(i)) = 0;
break;
end
end
end
end
else
accepted(ranking(i)) = 0;
end
end
if size(comp1,2) == 2
A0=convex_hull(1,end-1); 
B0=convex_hull(2,end-1); 
end
x = zeros(1, sum(accepted(:)));
y = zeros(1, sum(accepted(:)));
x_= zeros(1, sum(newInd(:)));
y_= zeros(1, sum(newInd(:)));
a1 = 0;
a1_new= 0;
count = 0; 
fp1 = fopen('extended_convex_hull', 'w');
[nothing, nothing] = unix(['cat /dev/null > extended_convex_hull_POSCARS']);  
fprintf(fp1,'# It contains the all the information in extendedConvexHull.pdf\n');
fprintf(fp1,'# X axis: X = y/(x+y) \n');
fprintf(fp1,'# Y axis: Y = (E(AxBy)-x*E(A)-y*E(B))/(x+y)\n');
fprintf(fp1,'# Fitness: its distance to the convex hull (eV/block)\n');
fprintf(fp1,'  ID   Compositions    Enthalpies     Volumes     Fitness   SYMM    X        Y     \n');
if size(comp,2) == 2  
fprintf(fp1,'                       (eV/atom)    (A^3/atom)   (eV/block)              (eV/atom) \n');
else
fprintf(fp1,'                       (eV/atom)    (A^3/atom)   (eV/block)              (eV/block) \n');
end
for i = 1 : Na
j = ranking(i);
if accepted(j)>0
if size(comp1,2) == 2
a1 = a1 + 1;
if size(numIons,2)==2
x(a1) = comp(j,2)/(comp(j,1) + comp(j,2));
y(a1) = (enth(j)*sum(comp(j,:)) - A0*sum(numIons(1,:))*comp1(j,1) - B0*sum(numIons(2,:))*comp1(j,2)) / sum(comp(j,:)); 
if newInd(j)==1
a1_new=a1_new+1;
x_(a1_new)=x(a1);
y_(a1_new)=y(a1);
end
else
x(a1) = comp1(j,2)/(comp1(j,1) + comp1(j,2));
y(a1) = (enth(j)*sum(comp(j,:)) - A0*sum(numIons(1,:))*comp1(j,1) - B0*sum(numIons(2,:))*comp1(j,2)) / sum(comp1(j,:)); 
if newInd(j)==1
a1_new=a1_new+1;
x_(a1_new)=x(a1);
y_(a1_new)=y(a1);
end
end
end
symg      = USPEX_STRUC.POPULATION(j).symg;
num       = USPEX_STRUC.POPULATION(j).numIons;
lat       = USPEX_STRUC.POPULATION(j).LATTICE;
coords    = USPEX_STRUC.POPULATION(j).COORDINATES;
order     = USPEX_STRUC.POPULATION(j).order;
numBlocks = USPEX_STRUC.POPULATION(j).numBlocks;
volume    = USPEX_STRUC.POPULATION(j).Vol/sum(num);
composition = sprintf('%3d',num);
shift=[4, 2, 1]; 
if size(composition,2)<11
composition=[composition,blanks(shift(length(num)))];
end
if size(comp1,2) == 2
fprintf(fp1,'%4d  [%11s]   %9.4f    %9.4f   %9.4f   %4d  %6.3f  %7.4f\n', ...
j, composition, enth(j), volume, fitness(j), symg, x(a1), y(a1));
else 
fprintf(fp1,'%4d  [%11s]   %9.4f    %9.4f   %9.4f   %4d\n', ...
j, composition, enth(j), volume, fitness(j), symg);
end
Write_POSCAR(atomType, j, symg, num, lat, coords);
[nothing, nothing] = unix([' cat POSCAR      >> extended_convex_hull_POSCARS']);
if count <  popSize
ratio = num/sum(num);
store = 0;
newcomposition=1;
for ind=1:size(POOL_STRUC.Composition_ratio,1)
if sum(abs(ratio-POOL_STRUC.Composition_ratio(ind,:))) < 0.01  
if POOL_STRUC.Composition_surviving(ind)<stopCrit 
E_ranking=POOL_STRUC.Composition_ranking(ind) + 1;  
if E_ranking == 1
enthalpy1 = enth(j);
POOL_STRUC.Composition_fitness(ind) = fitness(j);
if enthalpy1 < POOL_STRUC.Composition_Bestenthalpy(ind)
POOL_STRUC.Composition_Bestenthalpy(ind) = enthalpy1;   
POOL_STRUC.Composition_surviving(ind) = 1;
else
POOL_STRUC.Composition_surviving(ind) = POOL_STRUC.Composition_surviving(ind)+1;
end
end
if E_ranking < bestHM
store = 1;
POOL_STRUC.Composition_ranking(ind) = E_ranking;
end
else
POOL_STRUC.Composition_ranking(ind) = -1 ; 
end
newcomposition = 0;
break;
end
end
if newcomposition==1;  
POOL_STRUC.Composition_ratio(end+1,:)= ratio;
POOL_STRUC.Composition_ranking(end+1) = 1;
POOL_STRUC.Composition_surviving(end+1)= 1;
POOL_STRUC.Composition_fitness(end+1)= fitness(j);
POOL_STRUC.Composition_Bestenthalpy(end+1)= enth(j);
store = 1;
end
if store == 1  
count     = count + 1;
POOL_STRUC.POPULATION(count).LATTICE = lat;
POOL_STRUC.POPULATION(count).COORDINATES = coords;
POOL_STRUC.POPULATION(count).numIons = num;
POOL_STRUC.POPULATION(count).numBlocks = numBlocks;
POOL_STRUC.POPULATION(count).order = order;
POOL_STRUC.POPULATION(count).enthalpy = enth(j);
POOL_STRUC.POPULATION(count).Number = j;
end
end 
end
end
fclose(fp1);
makeFigure_convexhull_301(generation, convex_hull, comp1, [x; y], [x_; y_]);
function [fitness] = final_convex_hull(convex_hull, enthalpies, blocks, numIons)
N_P = size(blocks, 1);  
N_T = size(blocks, 2);  
N_atomType = size(numIons,2);
fitness = zeros(N_P,1);
convex_hull0 = Transform_Convex_Hull(convex_hull, numIons, 1);  
convex_hull_block = convex_hull0(:,1:N_T);
convex_hull_E     = convex_hull0(:,N_T+1);
for i = 1 : N_P
N_atom  = sum(blocks(i,:)*numIons);
N_Block = sum(blocks(i,:));
E = enthalpies(i)*N_atom/N_Block;
[tmp, samecomp, fitness(i)] = CheckDecomposition(convex_hull_block, ...
convex_hull_E, blocks(i,:), E);
end
fitness = round(fitness*10000)/10000;
