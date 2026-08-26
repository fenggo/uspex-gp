function createORG_Onthefly(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
hardCore = python_uspex(getPy, ['-f ' inputFile ' -b IonDistances -e EndDistances']);
if ~isempty(hardCore)
hardCore = str2num(hardCore);
if size(hardCore,1) == size(hardCore,2)   
for i = 1 : length(ORG_STRUC.atomType)
for j = i : length(ORG_STRUC.atomType)
ORG_STRUC.minDistMatrice(i,j) = hardCore(i,j);
ORG_STRUC.minDistMatrice(j,i) = hardCore(i,j);
end
end
end
end
numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if ~isempty(numParallelCalcs)
ORG_STRUC.numParallelCalcs = str2num(numParallelCalcs);
end
try
CenterminDistMatrice=python_uspex(getPy, '-f INPUT.txt -b MolCenters -e  EndMol');
if ~isempty(CenterminDistMatrice)
ORG_STRUC.CenterminDistMatrice = str2num(CenterminDistMatrice);
end
catch
end
