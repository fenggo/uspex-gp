function update_USPEX_INDIVIDUAL(POPULATION, resFolder, gen, atomType, flag)
global USPEX_STRUC
ID = POPULATION.Number;
if flag == 1
USPEX_STRUC.POPULATION(ID).gen         =  gen; 
USPEX_STRUC.POPULATION(ID).ToCount     =  0;   
elseif flag == 2
USPEX_STRUC.POPULATION(ID).ToCount     =  1;   
end
Common_var = {'COORDINATES','LATTICE','numIons','FINGERPRINT','order',...
'PressureTensor', 'Parents', 'symg','K_POINTS','Enthalpies', 'superCell',...
'coords0', 'lat0', 'PressureTensor0'};
for i = 1:length(Common_var)
eval(['USPEX_STRUC.POPULATION(ID).' Common_var{i} ' = POPULATION.' Common_var{i} ';']);
end
safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
