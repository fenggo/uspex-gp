function META_CreateORGStruc(inputFile)
global ORG_STRUC
ORG_STRUC = struct('dimension',{},'varcomp',{},'molecule',{},...       
'minDistMatrice',{}, 'numIons',{},...   
'atomType', {}, 'lattice', {}, 'coordinates', {}, ...              
'valences', {}, 'NvalElectrons', {}, 'goodBonds',{}, 'checkConnectivity', {},...
'numGenerations', {}, 'populationSize', {}, ...
'GaussianWidth', {}, 'GaussianHeight', {},...                      
'specificFolder',{}, 'homePath',{}, 'resFolder',{}, ...  
'pickUpYN', {}, 'pickUpGen', {}, 'pickUpFolder',{}, ...            
'RmaxFing', {}, 'deltaFing', {}, 'sigmaFing', {}, 'doFing', {}, ... 
'toleranceFing', {}, 'erf_table', {}, 'weight', {}, 'log_file', {},...         
'commandExecutable', {}, 'Kresol', {}, 'abinitioCode', {}, ...      
'platform', {}, 'remoteFolder',{}, 'SGtolerance',{}, ...
'maxErrors', {}, 'numProcessors', {},...
'numParallelCalcs',{}, 'doSpaceGroup',{},'ExternalPressure',{}, ...
'Softmode_Fre', {}, 'Softmode_num', {},'Softmode_eignFre', {},...   
'Softmode_eignVec', {}, 'Softmode_eignSupercells',{}, ...           
'FullRelax', {}, 'coorMutationDegree',{},'conv_till', {}, ...       
'maxAt',{}, 'correctionAngle',{},  ...  
'basicStructureNumber', {}, 'useBasicCell', {},...                  
'maxVectorLength', {}, 'repeatForStatistics', {},...                                           
'maxIncrease', {}, 'minVectorLength', {}, 'correctionLength', {});  
[nothing,homePath] = unix('pwd');
homePath(end) = [];
ORG_STRUC(1).homePath = homePath;
ORG_STRUC(1).USPEXPath= homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
if findstr(uspexmode,'exe')
[nothing, pathtemp] = unix('echo -e $USPEXPATH');;
ORG_STRUC(1).USPEXPath=pathtemp(1:end-1);
end
ORG_STRUC.dimension = 3;
ORG_STRUC.molecule  = 0;
ORG_STRUC.varcomp   = 1;
ORG_STRUC.specificFolder = 'Specific';
ORG_STRUC.log_file = 'OUTPUT.txt';
ORG_STRUC.basicStructureNumber = 0;
ORG_STRUC.checkConnectivity = 1; 
ORG_STRUC.maxErrors = 2;
ORG_STRUC.useBasicCell = 5;
ORG_STRUC.correctionAngle = 0;
ORG_STRUC.correctionLength = 0;
ORG_STRUC.repeatForStatistics = 1;
ORG_STRUC.erf_table = zeros(803,1);
for i = 1 : 803
ORG_STRUC.erf_table(i) = erf((i-402)/100);
end
ORG_STRUC.RmaxFing = 10;
ORG_STRUC.deltaFing = 0.08;
ORG_STRUC.sigmaFing = 0.03;
ORG_STRUC.toleranceFing = 0.008;
ORG_STRUC.toleranceBestHM = 0.02;
ORG_STRUC.Softmode_num = 0;
createORG_System(inputFile);
createORG_AbinitCode(inputFile);
createORG_Symmetry(inputFile);
createORG_Fingerprint(inputFile);
createORG_META(inputFile);
existfold = 1;
folderNum = 0;
while ~isempty(existfold)
folderNum = folderNum + 1;
ORG_STRUC.resFolder = ['results' num2str(folderNum)];
[nothing1,existfold] = mkdir(ORG_STRUC.resFolder);
end
[nothing,homePath] = unix('pwd');
homePath(end) = [];
ORG_STRUC.homePath = homePath;
writeORG();
