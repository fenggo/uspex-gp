function update_STUFF_old(inputFile, goodFrac, ranking)
global ORG_STRUC
if size(ORG_STRUC.numIons,2)==1
ORG_STRUC.fracPerm = 0;
ORG_STRUC.fracGene = ORG_STRUC.fracGene+ORG_STRUC.fracPerm;
end
howManyOffsprings     = round(ORG_STRUC.populationSize*ORG_STRUC.fracGene);
howManyPermutations   = round(ORG_STRUC.populationSize*ORG_STRUC.fracPerm);
howManyAtomMutations  = round(ORG_STRUC.populationSize*ORG_STRUC.fracAtomsMut);
howManyRand           = round(ORG_STRUC.populationSize*ORG_STRUC.fracRand);
howManyRotations      = round(ORG_STRUC.populationSize*ORG_STRUC.fracRotMut);      
howManyTransmutations = round(ORG_STRUC.populationSize*ORG_STRUC.fracTrans);       
howManySecSwitch      = round(ORG_STRUC.populationSize*ORG_STRUC.fracSecSwitch);   
howManyShiftBorder    = round(ORG_STRUC.populationSize*ORG_STRUC.fracShiftBorder); 
sum_offsprings = howManyOffsprings     + howManyPermutations  + ...
howManyTransmutations + howManyAtomMutations + ...
howManyRand           + howManyRotations     + ...
howManySecSwitch      + howManyShiftBorder;
if sum_offsprings > ORG_STRUC.populationSize
howManyOffsprings     = round((howManyOffsprings    *ORG_STRUC.populationSize)/sum_offsprings);
howManyPermutations   = round((howManyPermutations  *ORG_STRUC.populationSize)/sum_offsprings);
howManyAtomMutations  = round((howManyAtomMutations *ORG_STRUC.populationSize)/sum_offsprings);
howManyRand           = round((howManyRand          *ORG_STRUC.populationSize)/sum_offsprings);
howManyTransmutations = round((howManyTransmutations*ORG_STRUC.populationSize)/sum_offsprings);
howManyRotations      = round((howManyRotations     *ORG_STRUC.populationSize)/sum_offsprings);
howManySecSwitch      = round((howManySecSwitch     *ORG_STRUC.populationSize)/sum_offsprings);
howManyShiftBorder    = round((howManyShiftBorder   *ORG_STRUC.populationSize)/sum_offsprings);
sum_offsprings        = howManyOffsprings     + howManyPermutations  + ...
howManyTransmutations + howManyAtomMutations + ...
howManyRand           + howManyRotations     + ...
howManySecSwitch      + howManyShiftBorder;
end
howManyleft = max([ORG_STRUC.populationSize-sum_offsprings, 0]);
if ~ORG_STRUC.constLattice
howManyMutations = howManyleft;  
else
howManyOffsprings = howManyOffsprings + howManyleft;
howManyMutations = 0;  
end
ORG_STRUC.howManyMutations     = howManyMutations;
ORG_STRUC.howManyPermutations  = howManyPermutations;
ORG_STRUC.howManyAtomMutations = howManyAtomMutations;
ORG_STRUC.howManyRand          = howManyRand;
ORG_STRUC.howManyRotations     = howManyRotations;      
ORG_STRUC.howManyOffsprings    = howManyOffsprings;
ORG_STRUC.howManyTrans         = howManyTransmutations;
ORG_STRUC.howManySecSwitch     = howManySecSwitch;      
ORG_STRUC.howManyShiftBorder   = howManyShiftBorder;    
end
