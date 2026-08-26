function createORG_Molecule(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
checkMolecules = python_uspex(getPy, ['-f ' inputFile ' -b checkMolecules -c 1']);

if ~isempty(checkMolecules)
   ORG_STRUC.checkMolecules = str2num(checkMolecules);
   % disp(['checkMolecules ' checkMolecules ]);
end
read_molecules();
