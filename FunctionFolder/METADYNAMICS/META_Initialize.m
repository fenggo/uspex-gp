function META_Initialize()
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
[nothing,homePath] = unix('pwd');
homePath(end) = [];
USPEXPath= homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
if findstr(uspexmode,'exe')
[nothing,USPEXPath]= unix('echo -e $USPEXPATH');;
USPEXPath(end)=[];
end
META_Initialize_POP();
POP_STRUC.resFolder = ORG_STRUC.resFolder;
USPEX_STRUC = struct('POPULATION',{}, 'GENERATION',{}, 'atomType', {}, 'Fp_weight', {});
USPEX_STRUC(1).GENERATION = struct('Best_enth', {}, 'Best_enth_relaxed',{}, 'BestID',{});
USPEX_STRUC.GENERATION(1) = QuickStart(USPEX_STRUC.GENERATION);
USPEX_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons',{}, ...
'symg',{}, 'howCome',{},'order',{}, 'FINGERPRINT', {}, 'K_POINTS', {}, ...
'ToCount',{}, 'gen',{}, 'Enthalpies', {}, 'Parents',{},'Vol',{},...
'PressureTensor',{}, 'coords0', {}, 'lat0',{}, 'PressureTensor0',{});
USPEX_STRUC.POPULATION(1) = QuickStart(USPEX_STRUC.POPULATION);
USPEX_STRUC.atomType = ORG_STRUC.atomType;
USPEX_STRUC.Fp_weight = ORG_STRUC.weight;
safesave('Current_POP.mat', POP_STRUC);
safesave([POP_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
POP_STRUC.resFolder = ORG_STRUC.resFolder;
if ORG_STRUC.varcomp == 1
[nothing, nothing] = unix(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals']);
[nothing, nothing] = unix(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals']);
if ORG_STRUC.FullRelax > 0
[nothing, nothing] = unix(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals_relaxed']);
[nothing, nothing] = unix(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals_relaxed']);
end
else
[nothing, nothing] = unix(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals']);
[nothing, nothing] = unix(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals']);
if ORG_STRUC.FullRelax > 0
[nothing, nothing] = unix(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals_relaxed']);
[nothing, nothing] = unix(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals_relaxed']);
end
end
WriteGenerationStart();
