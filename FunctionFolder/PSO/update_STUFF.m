function update_STUFF(inputFile)
global ORG_STRUC
global POP_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if ~isempty(numParallelCalcs)
ORG_STRUC.numParallelCalcs = str2num(numParallelCalcs);
end
reoptOld = python_uspex(getPy, ['-f ' inputFile ' -b reoptOld -c 1']);
if ~isempty(reoptOld)
ORG_STRUC.reoptOld = str2num(reoptOld);
end
antiSeedsMax = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsMax -c 1']);
if isempty(antiSeedsMax)
antiSeedsMax = '0.000'; 
end
ORG_STRUC.antiSeedsMax = str2num(antiSeedsMax);
antiSeedsSigma = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsSigma -c 1']);
if isempty(antiSeedsSigma)
antiSeedsSigma = '0.001'; 
end
ORG_STRUC.antiSeedsSigma = str2num(antiSeedsSigma);
ordering = python_uspex(getPy, ['-f ' inputFile ' -b ordering_active -c 1']);
if isempty(ordering)
ordering = '1'; 
end
ORG_STRUC.ordering = str2num(ordering);
antiSeedsActivation = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsActivation -c 1']);
if isempty(antiSeedsActivation)
antiSeedsActivation = '1'; 
end
ORG_STRUC.antiSeedsActivation = str2num(antiSeedsActivation);
fracRand = python_uspex(getPy, ['-f ' inputFile ' -b fracRand -c 1']);
if isempty(fracRand)
fracRand = '0.1'; 
end
ORG_STRUC.fracRand = str2num(fracRand);
save ('Current_ORG.mat', 'ORG_STRUC')
