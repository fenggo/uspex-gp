function update_STUFF(inputFile, goodFrac, ranking)
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
createORG_Fingerprint(inputFile);
createORG_EA(inputFile);
createORG_Symmetry(inputFile);
createORG_Onthefly(inputFile);
N_pro = round((length(POP_STRUC.ranking))*ORG_STRUC.bestFrac);
if N_pro == 0
N_pro = 1;
end
howManyProliferate = min( [N_pro,length(POP_STRUC.ranking)-POP_STRUC.bad_rank ] );
if ~ORG_STRUC.AutoFrac
update_STUFF_old(inputFile, goodFrac, ranking);
else
Parent    = {'Random',    'Heredity', 'Permutate', 'TransMutate', ...
'LatMutate', 'softmutate', 'Rotate'};
fractions = {'fracRand', 'fracGene', 'fracPerm', 'fracTrans',...
'fracLatMut','fracAtomsMut','fracRotMut'};
howMany   = {'howManyRand', 'howManyOffsprings', 'howManyPermutations', 'howManyTrans',...
'howManyMutations','howManyAtomMutations','howManyRotations'};
numOperation = length(Parent); 
N = zeros(3,numOperation); 
f = zeros(1,numOperation); 
gen = POP_STRUC.generation;
if gen > 1
prev1_gen  = USPEX_STRUC.GENERATION(gen-1).ID;  
if gen == 2
prev2_gen  = 0;
else
prev2_gen  = USPEX_STRUC.GENERATION(gen-2).ID;  
end
for i = 1:howManyProliferate 
C_ID = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).Number;
accept = 1;
if C_ID > 0
for j = prev2_gen + 1 : prev1_gen
if ORG_STRUC.dimension == 3
ans = SameStructure_order(j, C_ID, USPEX_STRUC);  
else
ans = SameStructure(j, C_ID, USPEX_STRUC);
end
if ans
accept = 0;
break;
end
end
end
if accept
tmp = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).howCome;
for count = 1:numOperation
if ~isempty(findstr(tmp, Parent{count}))
N(1,count) = N(1,count) + 1;
break;
end
end
end
end
for i = 1:length(POP_STRUC.POPULATION) 
tmp = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).howCome;
for count = 1:numOperation
if ~isempty(findstr(tmp, Parent{count}))
N(2,count) = N(2,count) + 1;
break;
end
end
end
count = 0;
sum1 = 0;
for i = 1:numOperation
if N(2,i) > 0
N(1,i) = max(1, N(1,i)); 
X(i) = N(1,i)/N(2,i);
count = count + 1;
sum1 = sum1 + X(i);
end
end
X_mean = sum1/count;
for i = 1:numOperation
if N(2,i) >0
N(3,i) = round(N(1,i)*(X(i)/X_mean+1)/2);
end
end
f = N(3,:)/sum(N(3,:));
end
for i = 1:numOperation
eval([ 'f(1,i) = 0.55 * ORG_STRUC.' fractions{i} '+ 0.45 * f(1,i);' ]);
end
f(1,1) = max(f(1,1), 0.10);  
f(1,2) = max(f(1,2), 0.10);  
f(1,6) = max(f(1,6), 0.10);  
if ORG_STRUC.fracTrans> 0
f(1,4) = max(f(1,4), 0.05);
end
f = f/sum(f);
pop = ORG_STRUC.populationSize;
for i = 1:numOperation
eval( ['ORG_STRUC.' fractions{i} '= f(1,i);'] );
eval( ['ORG_STRUC.' howMany{i} '= round(pop*f(1,i));'] );
end
end
ORG_STRUC.tournament = zeros(howManyProliferate,1);
ORG_STRUC.tournament(howManyProliferate) = 1;
for loop = 2:howManyProliferate
ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2) + loop^2;
end
if POP_STRUC.generation == 0
if howManyProliferate > ORG_STRUC.initialPopSize
if ORG_STRUC.initialPopSize==0
ORG_STRUC.initialPopSize=1;
end
ORG_STRUC.tournament = zeros(ORG_STRUC.initialPopSize,1);
ORG_STRUC.tournament(ORG_STRUC.initialPopSize) = 1;
for loop = 2:ORG_STRUC.initialPopSize
ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2)+loop^2;
end
end
end
if size(ORG_STRUC.softMutOnly,1)<ORG_STRUC.numGenerations+2
ORG_STRUC.softMutOnly = zeros(ORG_STRUC.numGenerations + 1,1);
end
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');
fprintf(fp, [alignLine('VARIATION OPERATORS') '\n']);
fprintf(fp,'The fittest %2d percent of the population used to produce next generation\n', ORG_STRUC.bestFrac*100);
fprintf(fp,'    fraction of generation produced by heredity        :     %4.2f\n', ORG_STRUC.fracGene);
fprintf(fp,'    fraction of generation produced by random          :     %4.2f\n', ORG_STRUC.fracRand);
fprintf(fp,'    fraction of generation produced by softmutation    :     %4.2f\n', ORG_STRUC.fracAtomsMut);
fprintf(fp,'    fraction of generation produced by permutation     :     %4.2f\n', ORG_STRUC.fracPerm);
fprintf(fp,'    fraction of generation produced by latmutation     :     %4.2f\n', ORG_STRUC.fracLatMut);
fprintf(fp,'    fraction of generation produced by rotmutation     :     %4.2f\n', ORG_STRUC.fracRotMut);
fprintf(fp,'    fraction of generation produced by transmutation   :     %4.2f\n', ORG_STRUC.fracTrans);
fclose(fp);
safesave ('Current_ORG.mat', ORG_STRUC)
