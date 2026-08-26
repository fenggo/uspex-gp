function loops_seed = pick_Seeds()
global POP_STRUC
global ORG_STRUC
warning off
cd Seeds
loops_seed = 0;
if exist([ 'POSCARS_', num2str(POP_STRUC.generation) ])
if exist('POSCARS')
[nothing, nothing] = unix([ 'cat POSCARS >> POSCARS_' num2str(POP_STRUC.generation) ]);
end
[nothing, nothing] = unix([ 'cp POSCARS_' num2str(POP_STRUC.generation) ' POSCARS' ]);
end
if exist('POSCARS')
[fid,message] = fopen('POSCARS');
while 1
try
tmp = fgetl(fid); 
scale_factor = fgetl(fid); 
optlattice = fscanf(fid,'%g',[3,3]); % lattice vectors
candidate_lat = optlattice'*str2num(scale_factor);
tmp = fgetl(fid);  
atomType = fgetl(fid); 
ntyp = fgetl(fid); 
ntyp = str2num(ntyp);
atomType1 = zeros(1,length(ntyp));
c1 = findstr(atomType, ' ');
c = sort(str2num(['0 ' num2str(c1)]));
c(end+1) = length(atomType) + 1;
ind1 = 1;
for i = 2 : length(c)
if c(i-1)+1 > c(i)-1
continue
end
tmp = atomType(c(i-1)+1 : c(i)-1);
if ~isempty(str2num(tmp)) 
atomType1(ind1) = str2num(tmp);
else
for j = 1 : 105
if strcmp(lower(tmp), lower(elementFullName(j))) | strcmp(lower(tmp), lower(megaDoof(j)))
atomType1(ind1) = j;
break;
end
end
end
ind1 = ind1 + 1;
end
loops_seed = loops_seed + 1;
for i=1:length(atomType1)
if atomType1(i)==0
disp(['Error, elemental types are inconsistent when reading Seeds' num2str(loops_seed)]);
disp('Program STOPS');
quit
end
end
numIons = zeros(1,length(ORG_STRUC.atomType));
for i=1:length(ORG_STRUC.atomType)
for j=1:length(ntyp)
if atomType1(j)==ORG_STRUC.atomType(i)
numIons(i) = ntyp(j);
end
end
end
if ORG_STRUC.molecule == 1
comp=zeros(size(ORG_STRUC.numIons,2), length(ORG_STRUC.atomType));
for i=1:size(ORG_STRUC.numIons,2)
for j=1:length(ORG_STRUC.STDMOL(i).types)
for k=1:length(ORG_STRUC.atomType)
if ORG_STRUC.STDMOL(i).types(j)==k
comp(i,k)=comp(i,k)+1;
end
end
end
end
end
if ORG_STRUC.varcomp==0
if ORG_STRUC.molecule==1
numIons1 = ORG_STRUC.numIons*comp;
if sum(numIons1 ~= numIons)>0
disp(['the number of atoms is inconsistent when reading Seeds-' num2str(loops_seed)]);
disp('Program STOPS');
quit
else
numMols = ORG_STRUC.numIons;
end
elseif numIons ~= ORG_STRUC.numIons
disp(['the number of atoms is inconsistent when reading Seeds-' num2str(loops_seed)]);
disp('Program STOPS');
quit
end
else
if ORG_STRUC.molecule==1
numMols = round(numIons/comp);
numBlocks1 = numMols/ORG_STRUC.numIons;
else
numBlocks1 = numIons/ORG_STRUC.numIons;
end
if abs( sum(numBlocks1-round(numBlocks1)) ) > 0.01
disp(['Impossible to make combination ' num2str(numBlocks1) ' when reading Seeds-' num2str(loops_seed)]);
disp('Program STOPS');
quit
end
numBlocks = round(numBlocks1);
end
natom = sum(ntyp); 
tmp = fgetl(fid);
candidate_pop = fscanf(fid,'%g',[3,natom]);
candidate_pop = candidate_pop';
tmp = fgetl(fid);
Add = length(POP_STRUC.POPULATION) + 1;
if POP_STRUC.generation == 1
ORG_STRUC.initialPopSize = ORG_STRUC.initialPopSize + 1;
end
info_parents = struct('parent',{}, 'enthalpy', {});
info_parents(1).parent = [];
POP_STRUC.POPULATION(Add).Parents = info_parents;
POP_STRUC.POPULATION(Add).LATTICE = candidate_lat;
POP_STRUC.POPULATION(Add).COORDINATES = candidate_pop;
POP_STRUC.POPULATION(Add).numIons = numIons;
POP_STRUC.POPULATION(Add).howCome = 'Seeds';
if ORG_STRUC.molecule==1
POP_STRUC.POPULATION(Add).numMols = numMols;
[type, MtypeLIST, numIons]=GetPOP_MOL(numMols);
POP_STRUC.POPULATION(Add).typesAList = type;
POP_STRUC.POPULATION(Add).MtypeLIST = MtypeLIST;
readMOL(Add, 0);
end
if ORG_STRUC.varcomp==1
POP_STRUC.POPULATION(Add).numBlocks = numBlocks;
end
disp(['seed number ' num2str(loops_seed) ' has been successfully added']);
catch
break
end
end  
status = fclose(fid);
disp('end of pickup Seeds')
[nothing, nothing] = unix([ 'echo Generation:' num2str(POP_STRUC.generation) '>> ../' ORG_STRUC.resFolder '/Seeds_history' ]);
[nothing, nothing] = unix([ 'cat POSCARS >> ../' ORG_STRUC.resFolder '/Seeds_history' ]);
[nothing, nothing] = unix([ 'mv POSCARS POSCARS_' num2str(POP_STRUC.generation) ]);
end
cd ..
