function NEB_checkInitialStuct()
global ORG_STRUC
global POP_STRUC
numImages =     ORG_STRUC.numImages;
msgVCNEB=0;
errVCNEB=0;
if ORG_STRUC.CalcType==1 && ORG_STRUC.numImages<3
disp(' ');
disp([ 'numImages = ' num2str(ORG_STRUC.numImages)]);
disp(' ');
disp('ERROR:  At least THREE Images should be set for VC-NEB calculation ! VCNEB Code quite with error ...');
disp(' ');
exit;
end
if ORG_STRUC.CalcType==1 && (ORG_STRUC.optMethodCIDI > 0)
if abs(ORG_STRUC.optMethodCIDI)==1 
if length(ORG_STRUC.pickupImages)>1
disp(' ');
disp([ 'pickupImages = ' num2str(ORG_STRUC.pickupImages)]);
disp(' ');
disp('ERROR:  Less than ONE Image is needed for Climbing/Downing Image method ! VCNEB Code quite with error ...');
disp(' ');
exit;
end
end
if abs(ORG_STRUC.optMethodCIDI)==2
if isempty(ORG_STRUC.pickupImages)
disp(' ');
disp([ 'pickupImages = ' num2str(ORG_STRUC.pickupImages)]);
disp(' ');
disp('ERROR:  More than ONE Image is needed for Multi-Climbing/Downing Image method ! VCNEB Code quite with error ...');
disp(' ');
exit;
end
end
if (~isempty(find(ORG_STRUC.pickupImages==1)) || ~isempty(find(ORG_STRUC.pickupImages==numImages))) && (abs(ORG_STRUC.optVarImage)>0) 
disp(' ');
disp([ 'pickupImages = ' num2str(ORG_STRUC.pickupImages)]);
disp(' ');
disp('ERROR:  Please donot assign the initial and final Image as the Climbing/Downing Image! VCNEB Code quite with error ...');
disp(' ');
exit;
end
if ~isempty(ORG_STRUC.pickupImages) && (ORG_STRUC.optVarImage==1)
disp(' ');
disp([ 'pickupImages = ' num2str(ORG_STRUC.pickupImages)]);
disp(' ');
disp('ERROR:  Please donot use Variable-Image VCNEB method when running Climbing/Downing Images calculations! VCNEB Code quite with error ...');
disp(' ');
exit;
end
end
return;
msgVCNEB.CalcType(1) = [' Calculation Type:  The VC-NEB calculation'];
msgVCNEB.CalcType(2) = [' Calculation Type:  Strucutre Optimization  '];
msgVCNEB.CalcType(3) = [' Calculation Type:  Strucutre Optimization with Temperature (not yet)'];
msgVCNEB.CalcType(4) = [' Calculation Type:  Rotation avoiding method Calculation'];
msgVCNEB.optRelaxType(1) = [' Relaxation Type: Only atomic postions are optimized '];
return
[nothing, numIons]   = unix(['./getFromTo numIons EndNumIons  '    inputFile]);
if isempty(numIons)
disp('Error : NO ATOM NUMBERS, the VCNEB code cannot run ! ');
else
ORG_STRUC.numIons = str2num(numIons);
end
[nothing, atomType]  = unix(['./getFromTo atomType EndAtomType '   inputFile]);
if isempty(atomType)
disp('Warnning : NO ATOM TYPE');
else
ORG_STRUC.atomType = str2num(atomType);
end
[nothing, SuperCell] = unix(['./getFromTo superCell EndSuperCell ' inputFile]);
if isempty(SuperCell)
else
ORG_STRUC.SuperCell= str2num(SuperCell);
end
