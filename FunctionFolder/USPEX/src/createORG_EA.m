function createORG_EA(inputFile)
warning off
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
numGenerations = python_uspex(getPy, ['-f ' inputFile ' -b numGenerations -c 1']);
if ~isempty(numGenerations)
ORG_STRUC.numGenerations = str2num(numGenerations);
end
stopCrit = python_uspex(getPy, ['-f ' inputFile ' -b stopCrit -c 1']);
if isempty(stopCrit)
if ORG_STRUC.varcomp
ORG_STRUC.stopCrit = ORG_STRUC.maxAt; 
else
ORG_STRUC.stopCrit = sum(ORG_STRUC.numIons); 
end
if ORG_STRUC.molecule
ORG_STRUC.stopCrit = ORG_STRUC.stopCrit*3; 
end
else
ORG_STRUC.stopCrit = str2num(stopCrit);
end
populationSize = python_uspex(getPy, ['-f ' inputFile ' -b populationSize -c 1']);
if isempty(populationSize)
if ORG_STRUC.varcomp==1
pS = 10*round(2*ORG_STRUC.maxAt/10); 
else
pS = 10*round(2*sum(ORG_STRUC.numIons)/10); 
end
if pS < 10
pS = 10;
elseif pS > 60
pS = 60;
end
populationSize = num2str(pS);
end
ORG_STRUC.populationSize = str2num(populationSize);
bestFrac = python_uspex(getPy, ['-f ' inputFile ' -b bestFrac -c 1']);
if ~isempty(bestFrac)
ORG_STRUC.bestFrac = str2num(bestFrac);
end
AutoFrac = python_uspex(getPy, ['-f ' inputFile ' -b AutoFrac -c 1']);
if ~isempty(AutoFrac)
ORG_STRUC.AutoFrac = str2num(AutoFrac);
end
fracGene = python_uspex(getPy, ['-f ' inputFile ' -b fracGene -c 1']);
if ~isempty(fracGene)
ORG_STRUC.fracGene = str2num(fracGene);
end
fracRand = python_uspex(getPy, ['-f ' inputFile ' -b fracRand -c 1']);
if ~isempty(fracRand)
ORG_STRUC.fracRand = str2num(fracRand);
end
fracPerm = python_uspex(getPy, ['-f ' inputFile ' -b fracPerm -c 1']);
if isempty(fracPerm)
if (length(ORG_STRUC.atomType) > 1) & (ORG_STRUC.varcomp==0)
fracPerm = '0.1'; 
else
fracPerm = '0.0'; 
end
end
ORG_STRUC.fracPerm = str2num(fracPerm);
fracLatMut = python_uspex(getPy, ['-f ' inputFile ' -b fracLatMut -c 1']);
if isempty(fracLatMut)
if ORG_STRUC.constLattice == 0
fracLatMut = '0.1'; 
else
fracLatMut = '0.0'; 
end
end
ORG_STRUC.fracLatMut = str2num(fracLatMut);
fracAtomsMut = python_uspex(getPy, ['-f ' inputFile ' -b fracAtomsMut -c 1']);
if ~isempty(fracAtomsMut)
ORG_STRUC.fracAtomsMut = str2num(fracAtomsMut);
end
if ORG_STRUC.molecule==1 | ((ORG_STRUC.dimension==-4)  &&  (ORG_STRUC.molecule==0)  &&  (ORG_STRUC.varcomp==0))
fracRotMut = python_uspex(getPy, ['-f ' inputFile ' -b fracRotMut -c 1']);
if isempty(fracRotMut)
fracRotMut = '0.1'; 
end
ORG_STRUC.fracRotMut = str2num(fracRotMut);
else
ORG_STRUC.fracRotMut = 0;
end
if ORG_STRUC.dimension == -4
fracSecSwitch = python_uspex(getPy, ['-f ' inputFile ' -b fracSecSwitch -c 1']);
if isempty(fracSecSwitch)
fracSecSwitch = '0.0'; 
end
ORG_STRUC.fracSecSwitch = str2num(fracSecSwitch);
fracShiftBorder = python_uspex(getPy, ['-f ' inputFile ' -b fracShiftBorder -c 1']);
if isempty(fracShiftBorder)
fracShiftBorder = '0.0'; 
end
ORG_STRUC.fracShiftBorder = str2num(fracShiftBorder);
end
if ORG_STRUC.varcomp == 1 | ORG_STRUC.dimension == 2
fracTrans = python_uspex(getPy, ['-f ' inputFile ' -b fracTrans -c 1']);
if isempty(fracTrans)
if size(ORG_STRUC.numIons, 1) >1 | ORG_STRUC.dimension==2
fracTrans = '0.1'; 
else
fracTrans = '0'; 
end
end
ORG_STRUC.fracTrans = str2num(fracTrans);
else
ORG_STRUC.fracTrans = 0;
end
howManyMut = python_uspex(getPy, ['-f ' inputFile ' -b mutationDegree -c 1']);
if isempty(howManyMut)
howManyMut = 0;
for i = 1 : length(ORG_STRUC.atomType)
howManyMut = howManyMut + str2num(covalentRadius(ceil(ORG_STRUC.atomType(i))));
end
ORG_STRUC.howManyMut = 3*howManyMut/length(ORG_STRUC.atomType);
else
ORG_STRUC.howManyMut = str2num(howManyMut);
end
mutationRate = python_uspex(getPy, ['-f ' inputFile ' -b mutationRate -c 1']);
if ~isempty(mutationRate)
ORG_STRUC.mutationRate = str2num(mutationRate);
end
softMutOnly = python_uspex(getPy, ['-f ' inputFile ' -b softMutOnly -e EndSoftOnly']);
ORG_STRUC.softMutOnly = zeros(ORG_STRUC.numGenerations + 1,1);
c1 = findstr(softMutOnly, ' ');
c2 = findstr(softMutOnly, '-');
c = sort(str2num(['0 ' num2str(c1) ' ' num2str(c2)]));
c(end+1) = length(softMutOnly) + 1;
ind1 = 1;
for i = 2 : length(c)
if c(i-1)+1 > c(i)-1
continue
end
ind2 = str2num(softMutOnly(c(i-1)+1 : c(i)-1));
if ind2 == 0
ind1 = 1;
continue
end
if ~isempty(find(c2 == c(i-1)))
for j = ind1 : ind2
ORG_STRUC.softMutOnly(j) = 1;
end
else
ORG_STRUC.softMutOnly(ind2) = 1;
end
ind1 = ind2;
end
softMutTill = python_uspex(getPy, ['-f ' inputFile ' -b softMutTill -c 1']);
if isempty(softMutTill)
softMutTill = num2str(ORG_STRUC.numGenerations); % default
end
ORG_STRUC.softMutTill = str2num(softMutTill);
if ORG_STRUC.varcomp==1 | ORG_STRUC.dimension==2
howManyTrans = python_uspex(getPy, ['-f ' inputFile ' -b howManyTrans -c 1']);
if isempty(howManyTrans)
if size(ORG_STRUC.numIons, 1) >1 | ORG_STRUC.dimension==2
howManyTrans = '0.2'; 
else
howManyTrans = '0.0'; 
end
end
ORG_STRUC.howManyTrans = str2num(howManyTrans);
end
try
specificSwaps = python_uspex(getPy, ['-f ' inputFile ' -b specificSwaps -e EndSpecific']);
ORG_STRUC.specificSwaps = str2num(specificSwaps);
if isempty(ORG_STRUC.specificSwaps)
ORG_STRUC.specificSwaps = 0;
end
catch
ORG_STRUC.specificSwaps = 0;
end
howManySwaps = python_uspex(getPy, ['-f ' inputFile ' -b howManySwaps -c 1']);
if isempty(howManySwaps) 
if ORG_STRUC.varcomp | (ORG_STRUC.fracPerm == 0)
howManySwaps = '0';
elseif (size(ORG_STRUC.specificSwaps,1) == 1) & (size(ORG_STRUC.specificSwaps,2) == 1) 
howManySwaps = num2str(round(min(ORG_STRUC.numIons)/2));
else
howMS = 0;
for i = 1 : size(ORG_STRUC.specificSwaps, 1)
howMS = howMS + min(ORG_STRUC.numIons(ORG_STRUC.specificSwaps(i,1)), ORG_STRUC.numIons(ORG_STRUC.specificSwaps(i,2)))/2;
end
howManySwaps = num2str(round(howMS));
end
end
ORG_STRUC.howManySwaps = str2num(howManySwaps);
try
specificTrans = python_uspex(getPy, ['-f ' inputFile ' -b specificTrans -e EndTransSpecific']);
ORG_STRUC.specificTrans = str2num(specificTrans);
if isempty(ORG_STRUC.specificTrans)
ORG_STRUC.specificTrans = 0;
end
catch
ORG_STRUC.specificTrans = 0;
end
manyParents = python_uspex(getPy, ['-f ' inputFile ' -b manyParents -c 1']);
if ~isempty(manyParents)
ORG_STRUC.manyParents = str2num(manyParents);
end
numparents = python_uspex(getPy, ['-f ' inputFile ' -b numberparents -c 1']);
if ~isempty(numparents)
ORG_STRUC.numParents = str2num(numparents);
end
if ORG_STRUC.numParents > 2
ORG_STRUC.manyParents = 1;
end
minSlice = python_uspex(getPy, ['-f ' inputFile ' -b minSlice -c 1']);
if ~isempty(minSlice)
ORG_STRUC.minSlice = str2num(minSlice);
end
maxSlice = python_uspex(getPy, ['-f ' inputFile ' -b maxSlice -c 1']);
if ~isempty(maxSlice)
ORG_STRUC.maxSlice = str2num(maxSlice);
end
percSliceShift = python_uspex(getPy, ['-f ' inputFile ' -b percSliceShift -c 1']);
if ~isempty(percSliceShift)
ORG_STRUC.percSliceShift = str2num(percSliceShift);
end
maxDistHeredity = python_uspex(getPy, ['-f ' inputFile ' -b maxDistHeredity -c 1']);
if ~isempty(maxDistHeredity)
ORG_STRUC.maxDistHeredity = str2num(maxDistHeredity);
end
keepBestHM =  python_uspex(getPy, ['-f ' inputFile ' -b keepBestHM -c 1']);
if isempty(keepBestHM)
keepBestHM = num2str(round(0.15*ORG_STRUC.populationSize)); % default
end
ORG_STRUC.keepBestHM = str2num(keepBestHM);
if ORG_STRUC.keepBestHM < 1
ORG_STRUC.keepBestHM = 1;
end
dynamicalBestHM = python_uspex(getPy, ['-f ' inputFile ' -b dynamicalBestHM -c 1']);
if ~isempty(dynamicalBestHM)
ORG_STRUC.dynamicalBestHM = str2num(dynamicalBestHM);
end
reoptOld = python_uspex(getPy, ['-f ' inputFile ' -b reoptOld -c 1']);
if ~isempty(reoptOld)
ORG_STRUC.reoptOld = str2num(reoptOld);
end
