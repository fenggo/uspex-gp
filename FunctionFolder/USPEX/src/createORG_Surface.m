function createORG_Surface(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
reconstruct = python_uspex(getPy, ['-f ' inputFile ' -b reconstruct -c 1']);
if ~isempty(reconstruct)
ORG_STRUC.reconstruct = str2num(reconstruct);
end
bulk_stoi = python_uspex(getPy, ['-f ' inputFile ' -b StoichiometryStart -e StoichiometryEnd']);
if isempty(bulk_stoi)
ORG_STRUC.bulk_stoi = [1 1];
else
ORG_STRUC.bulk_stoi = str2num(bulk_stoi);
end
E_AB = python_uspex(getPy, ['-f ' inputFile ' -b E_AB -c 1']);
if ~isempty(E_AB)
ORG_STRUC.E_AB = str2num(E_AB);
end
E_A = python_uspex(getPy, ['-f ' inputFile ' -b Mu_A -c 1']);
if ~isempty(E_A)
ORG_STRUC.E_A = str2num(E_A);
end
E_B = python_uspex(getPy, ['-f ' inputFile ' -b Mu_B -c 1']);
if ~isempty(E_B)
ORG_STRUC.E_B = str2num(E_B);
end
