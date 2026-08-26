function Initialize()
global ORG_STRUC
global POP_STRUC
global PSO_STRUC
[nothing,homePath] = unix('pwd');
homePath(end) = [];
USPEXPath= homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
if findstr(uspexmode,'exe')
[nothing,USPEXPath]= unix('echo -e $USPEXPATH');;
USPEXPath(end)=[];
end
if ~exist('Seeds','dir')
[nothing, nothing] = unix('mkdir Seeds');
end
if ~exist('AntiSeeds','dir')
[nothing, nothing] = unix('mkdir AntiSeeds');
end
Initialize_POP_STRUC_PSO();
POP_STRUC.resFolder = ORG_STRUC.resFolder;
PSO_STRUC = struct('POPULATION',{}, 'GENERATION',{}, 'atomType', {}, 'Fp_weight', {});
PSO_STRUC(1).GENERATION = struct('quasiEntropy', {}, 'convex_hull', {}, 'composEntropy',{}, 'BestID',{});
PSO_STRUC.GENERATION(1) = QuickStart(PSO_STRUC.GENERATION);
PSO_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons',{}, ...
'symg',{}, 'howCome',{},'order',{}, 'Fitness', {}, 'FINGERPRINT', {}, 'K_POINTS', {}, ...
'ToCount',{},'S_order',{}, 'gen',{}, 'struc_entr',{}, 'Enthalpies', {}, 'Parents',{},'Vol',{},...
'cell',{},'Surface_numIons',{},'numBlocks',{},'gap','hardness','mag_moment','dielectric_tensor');
PSO_STRUC.POPULATION(1) = QuickStart(PSO_STRUC.POPULATION);
USPEX_STRUC.atomType = ORG_STRUC.atomType;
USPEX_STRUC.Fp_weight = ORG_STRUC.weight;
safesave('Current_POP.mat', POP_STRUC);
safesave([POP_STRUC.resFolder '/PSO.mat'], PSO_STRUC);
[nothing, nothing] = unix(['echo "Gen   ID    Origin   Composition    Enthalpy   Volume  Density   Fitness   KPOINTS  SYMM  Q_entr A_order S_order" >>' ORG_STRUC.resFolder '/Individuals']);
[nothing, nothing] = unix(['echo "                                      (eV)     (A^3)   (g/cm^3)" >>' ORG_STRUC.resFolder '/Individuals']);
[nothing, nothing] = unix(['echo "  ID   Origin     Composition  Enthalpy(eV)  Volume(A^3)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
[nothing, nothing] = unix(['echo " ID    Origin    Enthalpy   Parent-E   Parent-ID" >>' ORG_STRUC.resFolder '/origin']);
WriteGenerationStart();
