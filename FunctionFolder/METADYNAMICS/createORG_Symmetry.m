function createORG_Symmetry(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
doSpaceGroup = python_uspex(getPy, ['-f ' inputFile ' -b doSpaceGroup -c 1']);
if isempty(doSpaceGroup)
if (ORG_STRUC.dimension==0) | (ORG_STRUC.dimension==2) | (ORG_STRUC.dimension==-3)
doSpaceGroup = '0'; 
else
doSpaceGroup = '1'; 
end
end
ORG_STRUC.doSpaceGroup = str2num(doSpaceGroup);
symmetrize = python_uspex(getPy, ['-f ' inputFile ' -b symmetrize -c 1']);
if ~isempty(symmetrize)
ORG_STRUC.symmetrize = str2num(symmetrize);
end
SGtolerance = python_uspex(getPy, ['-f ' inputFile ' -b SymTolerance -c 1']);
if ~isempty(SGtolerance)
SGtolerance = deblank(SGtolerance);
if strcmp(lower(SGtolerance), 'high')
ORG_STRUC.SGtolerance = 0.05;
elseif strcmp(lower(SGtolerance), 'medium')
ORG_STRUC.SGtolerance = 0.10;
elseif strcmp(lower(SGtolerance), 'low')
ORG_STRUC.SGtolerance = 0.20;
else 
ORG_STRUC.SGtolerance = str2num(SGtolerance);
end
end
if isempty(ORG_STRUC.SGtolerance) 
ORG_STRUC.SGtolerance = 0.10;
end
sym_coef = python_uspex(getPy, ['-f ' inputFile ' -b constraint_enhancement -c 1']);
if ~isempty(sym_coef)
ORG_STRUC.sym_coef = str2num(sym_coef);
end
