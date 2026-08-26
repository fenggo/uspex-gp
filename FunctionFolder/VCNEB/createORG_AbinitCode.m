function createORG_AbinitCode(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
ORG_STRUC.abinitioCode = python_uspex(getPy, ['-f ' inputFile ' -b abinitioCode -e ENDabinit'], 1);
whichCluster = python_uspex(getPy, ['-f ' inputFile ' -b whichCluster -c 1']);
if ~isempty(whichCluster)
if ~isempty(str2num(whichCluster)) 
ORG_STRUC.platform=abs(str2num(whichCluster));
else
whichCluster(end) = [];
if strcmp(whichCluster, 'nonParallel')
ORG_STRUC.platform = 0;
elseif strcmp(whichCluster, 'CFN')
ORG_STRUC.platform = 3;
elseif strcmp(whichCluster, 'QSH')
ORG_STRUC.platform = 4;
elseif strcmp(whichCluster, 'QSH2')
ORG_STRUC.platform = 5;
elseif strcmp(whichCluster, 'xservDE')
ORG_STRUC.platform = 6;
elseif strcmp(whichCluster, 'MIPT')
ORG_STRUC.platform = 7;
end
end
end
if ORG_STRUC.platform > 7
disp(['you did not chose a valid platform. Program STOPS...']);
fprintf(fp,'you did not chose a valid platform. Program STOPS...\n');
quit;
end
Kresol = python_uspex(getPy, ['-f ' inputFile ' -b KresolStart -e Kresolend']);
if isempty(Kresol) 
Ltmp = length(ORG_STRUC.abinitioCode);
ORG_STRUC.Kresol = 0.2*ones(1, Ltmp);
if Ltmp > 1
for i = 1 : Ltmp
ORG_STRUC.Kresol(i) = 0.2 - (i-1)*0.12/(Ltmp-1);
end
end
else
ORG_STRUC.Kresol = str2num(Kresol);
end
remoteFolder1 = python_uspex(getPy, ['-f ' inputFile ' -b remoteFolder -c 1']);
ORG_STRUC.remoteFolder = remoteFolder1(1:length(remoteFolder1)-1);
numProcessors = python_uspex(getPy, ['-f ' inputFile ' -b numProcessors -e EndProcessors']);
if isempty(numProcessors)
ORG_STRUC.numProcessors = ones(1, length(ORG_STRUC.abinitioCode))*8; 
else
ORG_STRUC.numProcessors = str2num(numProcessors);
end
maxErrors = python_uspex(getPy, ['-f ' inputFile ' -b maxErrors -c 1']);
if ~isempty(maxErrors)
ORG_STRUC.maxErrors=str2num(maxErrors);
end
numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if ~isempty(numParallelCalcs)
ORG_STRUC.numParallelCalcs = str2num(numParallelCalcs);
end
if ORG_STRUC.platform==0
commandExecutable = python_uspex(getPy, ['-f ' inputFile ' -b commandExecutable -e EndExecutable']);
helper = double(commandExecutable);
tempIND = find (helper==10);
helpIND = zeros(length(tempIND)+1,1);
helpIND(2:end)=tempIND';
for listLoop = 1:length(helpIND)-1
listCommand{listLoop}=commandExecutable(helpIND(listLoop)+1:helpIND(listLoop+1)-1);
end
for commandLoop = 1: length(ORG_STRUC.abinitioCode)
if ORG_STRUC.abinitioCode(commandLoop) == 0
tempCom{commandLoop} = '0';     
elseif listLoop < length(ORG_STRUC.abinitioCode)
tempCom{commandLoop} = listCommand{1};
else
tempCom{commandLoop} = listCommand{commandLoop};
end
end
ORG_STRUC.commandExecutable=tempCom;
end
