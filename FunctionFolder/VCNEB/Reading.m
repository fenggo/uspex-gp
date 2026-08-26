function Error = Reading(code, Ind_No, indic)
global POP_STRUC
global ORG_STRUC
Error = 0;
cd ([ ORG_STRUC.homePath '/CalcFold' num2str(indic)])
Gen = POP_STRUC.generation;
Step = POP_STRUC.POPULATION(Ind_No).Step;
numIons = ORG_STRUC.numIons;
maxErrors = 3;
Const_Lat = 1;
ID = ['Gen' num2str(Gen) '-Ind' num2str(Ind_No) '-Step' num2str(Step)]; %For ERROR output
TotalStep = length([ORG_STRUC.abinitioCode]);
GoodBad = Read_AbinitCode(code, 0, ID);  
if GoodBad==1
if POP_STRUC.POPULATION(Ind_No).Error <= maxErrors;
POP_STRUC.POPULATION(Ind_No).Enthalpy(Step)   = Read_AbinitCode(code, 1, ID);
POP_STRUC.POPULATION(Ind_No).cellStressMatirx = Read_AbinitCode(code, 2, ID);
origReadingForce    = Read_AbinitCode(code, 7, ID);
POP_STRUC.POPULATION(Ind_No).atomForcesMatrix = origReadingForce;
if     code==1  
elseif code==3  
L=POP_STRUC.POPULATION(Ind_No).LATTICE;
POP_STRUC.POPULATION(Ind_No).atomForcesMatrix = -( (L'*L)*L^(-1)*origReadingForce' )'; 
elseif code==7  
elseif code==8  
end
POP_STRUC.POPULATION(Ind_No).Step = Step + 1;
POP_STRUC.POPULATION(Ind_No).Error=0;
end
Clean_AbinitCode(code);
else
Error = 1;
if (code~=1) & (code~=8) & (code~=9)  
POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
else
POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 1; 
end
disp(['PROBLEM_reading Structure ' num2str(Ind_No)]);
end
cd ..
