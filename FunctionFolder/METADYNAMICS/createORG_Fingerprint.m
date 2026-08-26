function createORG_Fingerprint(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
RmaxFing =  python_uspex(getPy, ['-f ' inputFile ' -b RmaxFing -c 1']);
if ~isempty(RmaxFing)
ORG_STRUC.RmaxFing = str2num(RmaxFing);
end
if ORG_STRUC.RmaxFing > 0
ORG_STRUC.doFing = 1;
else
ORG_STRUC.doFing = 0;
end
deltaFing = python_uspex(getPy, ['-f ' inputFile ' -b deltaFing -c 1']);
if ~isempty(deltaFing)
ORG_STRUC.deltaFing = str2num(deltaFing);
end
sigmaFing = python_uspex(getPy, ['-f ' inputFile ' -b sigmaFing -c 1']);
if ~isempty(sigmaFing)
ORG_STRUC.sigmaFing = str2num(sigmaFing);
end
toleranceFing = python_uspex(getPy, ['-f ' inputFile ' -b toleranceFing -c 1']);
if ~isempty(toleranceFing)
ORG_STRUC.toleranceFing = str2num(toleranceFing);
end
toleranceBestHM = python_uspex(getPy, ['-f ' inputFile ' -b toleranceBestHM -c 1']);
if ~isempty(toleranceBestHM)
ORG_STRUC.toleranceBestHM = str2num(toleranceBestHM);
end
antiSeedsMax = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsMax -c 1']);
if ~isempty(antiSeedsMax)
ORG_STRUC.antiSeedsMax = str2num(antiSeedsMax);
end
antiSeedsSigma = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsSigma -c 1']);
if ~isempty(antiSeedsSigma)
ORG_STRUC.antiSeedsSigma = str2num(antiSeedsSigma);
end
ordering = python_uspex(getPy, ['-f ' inputFile ' -b ordering_active -c 1']);
if ~isempty(ordering)
ORG_STRUC.ordering = str2num(ordering);
end
antiSeedsActivation = python_uspex(getPy, ['-f ' inputFile ' -b anti_seeds_activation -c 1']);
if ~isempty(antiSeedsActivation)
ORG_STRUC.antiSeedsActivation = str2num(antiSeedsActivation);
end
if ORG_STRUC.molecule == 1
numIons = ORG_STRUC.numMols;
else
numIons = ORG_STRUC.numIons;
end
L = length(numIons);
S = 0;
ORG_STRUC.weight = zeros(L*L,1);
for i = 1:L
for j = 1:L
ind = (i-1)*L+j;
ORG_STRUC.weight(ind) = (numIons(i)*numIons(j));
S = S + (numIons(i)*numIons(j));
end
end
ORG_STRUC.weight = ORG_STRUC.weight/S;
