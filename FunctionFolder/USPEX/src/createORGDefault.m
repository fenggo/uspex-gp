function createORGDefault()
global ORG_STRUC

[nothing,homePath] = unix('pwd');
homePath(end) = [];
ORG_STRUC(1).homePath = homePath;
ORG_STRUC(1).USPEXPath= homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');

if findstr(uspexmode,'exe')
   [nothing, pathtemp] = unix('echo -e $USPEXPATH');;
   ORG_STRUC(1).USPEXPath=pathtemp(1:end-1);
end

ORG_STRUC.fixRndSeed   = 0;
ORG_STRUC.collectForces= 0 ;
ORG_STRUC.PhaseDiagram = 0;
ORG_STRUC.checkMolecules = 1;
ORG_STRUC.symmetrize = 0;
ORG_STRUC.cor_dir = 1;
ORG_STRUC.correlation_coefficient = 1;
ORG_STRUC.averageFitness = 1000000;
ORG_STRUC.constLattice = 0;
ORG_STRUC.varcomp = 0;
ORG_STRUC.ldaU    = 0;
ORG_STRUC.softMutOnly = 0;
ORG_STRUC.splitN = 0;
ORG_STRUC.vacuumSize = 0;
ORG_STRUC.maxErrors = 2;
ORG_STRUC.sym_coef = 1;
ORG_STRUC.chargeNeutrality = 0;
ORG_STRUC.splitInto = 0;
ORG_STRUC.repeatForStatistics = 1;
ORG_STRUC.ordering = 1;
ORG_STRUC.minAngle = 55;   
ORG_STRUC.minDiagAngle = 30;  
ORG_STRUC.maxErrors = 2;
ORG_STRUC.bestFrac = 0.7;        
ORG_STRUC.alpha = -1; 
ORG_STRUC.mutationRate = 0.5; 
ORG_STRUC.numGenerations = 100; 
ORG_STRUC.percSliceShift = 1.0; 
ORG_STRUC.dynamicalBestHM = 2; 
ORG_STRUC.reoptOld = 0;  
ORG_STRUC.pickUpYN = 0;
ORG_STRUC.pickUpGen = 1;
ORG_STRUC.pickUpFolder = 1;
ORG_STRUC.pickUpNCount = 0;
ORG_STRUC.remote = 0;
ORG_STRUC.numParallelCalcs = 1;
ORG_STRUC.pluginType = 0;         
ORG_STRUC.startNextGen = 1;       
ORG_STRUC.currentGenDone = 0;     
ORG_STRUC.log_file = 'OUTPUT.txt';
ORG_STRUC.specificFolder = 'Specific';
ORG_STRUC.abinitioCode = 1;
ORG_STRUC.platform = 0;
ORG_STRUC.optType = 1;  
ORG_STRUC.opt_sign = 1; 
ORG_STRUC.dimension = 3;
ORG_STRUC.molecule  = 0;
ORG_STRUC.varcomp   = 0;
ORG_STRUC.maxDistHeredity = 0.5;
ORG_STRUC.symmetrize = 0;
ORG_STRUC.SGtolerance = 0.10;
ORG_STRUC.sym_coef = 1;
ORG_STRUC.RmaxFing = 10;
ORG_STRUC.deltaFing = 0.08;
ORG_STRUC.sigmaFing = 0.03;
ORG_STRUC.toleranceFing = 0.008;
ORG_STRUC.toleranceBestHM = 0.02;
ORG_STRUC.erf_table = zeros(803,1);

for i = 1 : 803
    ORG_STRUC.erf_table(i) = erf((i-402)/100);
end

ORG_STRUC.antiSeedsMax = 0;
ORG_STRUC.antiSeedsSigma = 0.001;
ORG_STRUC.antiSeedsActivation = 5000;
ORG_STRUC.reconstruct = 1;
ORG_STRUC.thicknessS  = 2.0;
ORG_STRUC.thicknessB  = 3.0;
ORG_STRUC.manyParents = 0;
ORG_STRUC.numParents = 2;
ORG_STRUC.minSlice = 2;
ORG_STRUC.maxSlice = 8;
ORG_STRUC.AutoFrac = 0;
% --- RL operator-bandit (L2) defaults ---
ORG_STRUC.RL_strategy    = 0;     % 0=off, 1=L1, 2=L2
ORG_STRUC.RL_forget      = 0.9;   % forgetting factor lambda
ORG_STRUC.RL_w1          = 1.0;   % reward weight: elite membership
ORG_STRUC.RL_w2          = 0.5;   % reward weight: enthalpy gain over parents
ORG_STRUC.RL_armCap      = 0.6;   % single-arm fraction cap
ORG_STRUC.RL_kappa       = 1.0;   % continuous-TS exploration scale
ORG_STRUC.RL_autoExport  = 1;     % export posterior at finish (default on)
ORG_STRUC.rotationAngleMax = pi/2;  % promoted from Rotation_310 hardcode
ORG_STRUC.translationMax   = 0.5;   % promoted from Rotation_310 hardcode
ORG_STRUC.softStepScale    = 1.0;   % softmode step-length multiplier
% --- end RL defaults ---
ORG_STRUC.fracGene = 0.5;
ORG_STRUC.fracRand = 0.2;
ORG_STRUC.fracAtomsMut = 0.1;
ORG_STRUC.fracLatMut = 0.1;
ORG_STRUC.fracRotMut = 0.1;
ORG_STRUC.fracTrans  = 0.1;
ORG_STRUC.fracPerm   = 0.1;
ORG_STRUC.mlPeerFilter = 0;
