function Reading(code, Ind_No, indic)
global POP_STRUC
global ORG_STRUC
cd (['CalcFold' num2str(indic)])
TotalStep = length([ORG_STRUC.abinitioCode]);
Step = POP_STRUC.POPULATION(Ind_No).Step;
Gen=POP_STRUC.generation;
ID = ['Gen' num2str(Gen) '-Ind' num2str(Ind_No) '-Step' num2str(Step)]; %For ERROR output
Const_Lat = 1; 
GoodBad = Read_AbinitCode(code, 0, ID);  
if GoodBad==1
[POP_STRUC.POPULATION(Ind_No).COORDINATES,POP_STRUC.POPULATION(Ind_No).LATTICE] = Read_Structure(code, Const_Lat);
POP_STRUC.POPULATION(Ind_No).Enthalpies(Step) = Read_AbinitCode(code, 1, ID);
POP_STRUC.POPULATION(Ind_No).PressureTensor   = Read_AbinitCode(code, 2, ID);
POP_STRUC.POPULATION(Ind_No).Step = Step + 1;
Clean_AbinitCode(code);
else
if (code~=1) & (code~=8) 
POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
else
POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 1; 
end
[nothing, nothing] = unix(['echo PROBLEM_reading Structure' num2str(Ind_No)]);
[nothing, nothing] = unix('pwd');
end
cd ..
