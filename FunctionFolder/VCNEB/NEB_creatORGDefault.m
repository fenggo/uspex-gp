function memStruct = NEB_creatORGDefault() 
global ORG_STRUC
memStruct.numeral = {
'numImages';
'CalcType';      'optRelaxType';   'optVarImage';      'optVarK';   
'optimizerType';  'optReadImages'; 'optMethodCIDI';    'optFreezing';  
'ionCh';         'whetherConstraint';
'dimension';                     
'Temperature';   'SuperCell'; 'ExternalPressure';
'dt';           'K_min';         'K_max';   'Kconstant';      
'VarPathLength'; 
'startCIDIStep';  
'ConvThreshold';     
'numSteps';   'FormatType';      'PrintStep';
'SGtolerance';
'continuedNow';
'pickUpYN';     'pickUpGen'; 'pickUpFolder'; 
'submitCount';  
'remoteRegime';              'numProcessors';     'numParallelCalcs';  
'UnitType';  'Ry2eV';        'KBar2GPa';  'fs';            'bohr2A';        'amu2eVs2_A2';   'THz2cm_1'; 
'read_Handles';  'write_Handles'; 
'platform';
'repeatForStatistics';
'maxErrors';
'RmaxFing'; 'deltaFing'; 'sigmaFing'; 'toleranceFing';
};
memStruct.array = {
'atomType';     'atomSymbol'; 
'numIons';      'numSpecies';
'SuperCell';
'whichCI';      'whichDI';     'pickupImages'; 
'Kresol'; 
'abinitioCode';
};
memStruct.string = {
'username';  
'whichCluster';   'wallTime';
'remoteFolder';   'specificFolder'; 
'inputFile';
'commandExecutable';
'resFolder';     'homePath'; 'USPEXPath';
'log_file';
};
ORG_STRUC = struct( 'dimension', 3 );
for  i = 1:length(memStruct.numeral)
ORG_STRUC = setfield(ORG_STRUC, memStruct.numeral{i}, 0 );
end
for  i = 1:length(memStruct.array)
ORG_STRUC = setfield(ORG_STRUC, memStruct.array{i}, '' );
end
for  i = 1:length(memStruct.string)
ORG_STRUC = setfield(ORG_STRUC, memStruct.string{i}, '' );
end
[nothing,homePath] = unix('pwd');
homePath(end) = [];
ORG_STRUC(1).homePath = homePath;
ORG_STRUC(1).USPEXPath= homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
if findstr(uspexmode,'exe')
[nothing, pathtemp] = unix('echo -e $USPEXPATH');
ORG_STRUC(1).USPEXPath=pathtemp(1:end-1);
end
ORG_STRUC.platform = 0;    
ORG_STRUC.CalcType = 1; 
ORG_STRUC.optRelaxType = 3;
ORG_STRUC.optimizerType = 1;
ORG_STRUC.optReadImages = 1;
ORG_STRUC.numImages = 1; 
ORG_STRUC.numSteps = 10; 
ORG_STRUC.ConvThreshold = 3e-3; 
ORG_STRUC.dt = 0.05;
ORG_STRUC.optVarK = 1;
ORG_STRUC.optVarImage = 1;
ORG_STRUC.VarPathLength = 0;
ORG_STRUC.K_min = 5;
ORG_STRUC.K_max = 5;
ORG_STRUC.Kconstant=5;
ORG_STRUC.optFreezing = 0;
ORG_STRUC.optMethodCIDI = 0;
ORG_STRUC.whichCI = 0;
ORG_STRUC.whichDI = 0;
ORG_STRUC.startCIDIStep = 100000; 
ORG_STRUC.ExternalPressure = 0;
ORG_STRUC.Temperature = 0;
ORG_STRUC.SuperCell = [1 1 1];
ORG_STRUC.UnitType = 1;
ORG_STRUC.PrintStep = 1; 
ORG_STRUC.abinitioCode = 1; 
ORG_STRUC.Kresol = 0.075;
ORG_STRUC.SGtolerance = 0.10;
ORG_STRUC.dimension  = 3;
ORG_STRUC.bohr2A     = 0.529177249;
ORG_STRUC.Ry2eV      = 13.6056981;
ORG_STRUC.KBar2GPa   = 10;
ORG_STRUC.fs         = 10^(-15);
ORG_STRUC.amu2eVs2_A2= 1.0364269*10^(-28);
ORG_STRUC.THz2cm_1   = (10^10)/29979245800 ;
ORG_STRUC.log_file = 'OUTPUT.txt';
ORG_STRUC.specificFolder = 'Specific';
ORG_STRUC.E_IniImage = 10E+9;
ORG_STRUC.fp_IniImage=[];
ORG_STRUC.E_FinImage = 10E+9;
ORG_STRUC.fp_FinImage=[];
ORG_STRUC.RmaxFing = 10;
ORG_STRUC.deltaFing = 0.04;
ORG_STRUC.sigmaFing = 0.03;
ORG_STRUC.toleranceFing = 0.005;
ORG_STRUC.erf_table = zeros(803,1);
for i = 1 : 803
ORG_STRUC.erf_table(i) = erf((i-402)/100);
end
ORG_STRUC.varcomp = 0;  
ORG_STRUC.molecule= 0; 
ORG_STRUC.submitCount = 0;
