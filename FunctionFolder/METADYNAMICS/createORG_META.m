function createORG_META(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
useBasicCell = python_uspex(getPy, ['-f ' inputFile ' -b retryBasicEvery -c 1']);
if ~isempty(useBasicCell)
ORG_STRUC.useBasicCell = str2num(useBasicCell); 
end
numGenerations = python_uspex(getPy, ['-f ' inputFile ' -b numGenerations -c 1']);
if isempty(numGenerations)
numGenerations = '100'; 
end
ORG_STRUC.numGenerations = str2num(numGenerations);
populationSize = python_uspex(getPy, ['-f ' inputFile ' -b populationSize -c 1']);
if isempty(populationSize)
pS = 10*round(2*sum(ORG_STRUC.numIons)/10); 
if pS < 10
pS = 10;
elseif pS > 60
pS = 60;
end
populationSize = num2str(pS);
end
ORG_STRUC.populationSize = str2num(populationSize);
maxIncrease = python_uspex(getPy, ['-f ' inputFile ' -b maxCellIncrease -c 1']);
if isempty(maxIncrease)
ORG_STRUC.maxIncrease = 5; 
else
ORG_STRUC.maxIncrease = str2num(maxIncrease); 
end
maxVectorLength = python_uspex(getPy, ['-f ' inputFile ' -b maxVectorLength -c 1']);
if isempty(maxVectorLength)
lat_tmp = latConverter(ORG_STRUC.lattice);
ORG_STRUC.maxVectorLength = 3*max(lat_tmp(1:3)); 
else
ORG_STRUC.maxVectorLength = str2num(maxVectorLength);
end
howManyMut = python_uspex(getPy, ['-f ' inputFile ' -b mutationDegree -c 1']);
if isempty(howManyMut)
howManyMut = 0;
for i = 1 : length(ORG_STRUC.atomType)
howManyMut = howManyMut + str2num(covalentRadius(ORG_STRUC.atomType(i)));
end
ORG_STRUC.howManyMut = 3*howManyMut/length(ORG_STRUC.atomType);
else
ORG_STRUC.howManyMut = str2num(howManyMut);
end
if ORG_STRUC.varcomp
maxAt = python_uspex(getPy, ['-f ' inputFile ' -b maxAt -c 1']);
if isempty(maxAt)
maxAt = '100'; 
end
ORG_STRUC.maxAt = str2num(maxAt);
else
ORG_STRUC.maxAt = sum(ORG_STRUC.numIons);
end
minVectorLength = python_uspex(getPy, ['-f ' inputFile ' -b minVectorLength -c 1']);
if isempty(minVectorLength)
minVectorLength = 0;
for i = 1 : length(ORG_STRUC.atomType)
if minVectorLength < 2*str2num(covalentRadius(ORG_STRUC.atomType(i)))
minVectorLength = 2*str2num(covalentRadius(ORG_STRUC.atomType(i)));
end
end
ORG_STRUC.minVectorLength = minVectorLength; 
else
ORG_STRUC.minVectorLength = str2num(minVectorLength);
end
GaussianWidth = python_uspex(getPy, ['-f ' inputFile ' -b GaussianWidth -c 1']);
if isempty(GaussianWidth)
lat1 = latConverter(lat);
GaussianWidth = 0.15*min(lat1(1:3));
disp(['Gaussian width = ' num2str(GaussianWidth)]);
disp(' ');
end
ORG_STRUC.GaussianWidth = str2num(GaussianWidth);
GaussianHeight = python_uspex(getPy, ['-f ' inputFile ' -b GaussianHeight -c 1']);
if isempty(GaussianHeight)
lat1 = latConverter(lat);
GaussianHeight = 1000*ORG_STRUC.GaussianWidth^2*min(lat1(1:3));
end
ORG_STRUC.GaussianHeight = str2num(GaussianHeight);
useBasicCell = python_uspex(getPy, ['-f ' inputFile ' -b retryBasicEvery -c 1']);
if ~isempty(useBasicCell)
ORG_STRUC.useBasicCell = str2num(useBasicCell); 
end
