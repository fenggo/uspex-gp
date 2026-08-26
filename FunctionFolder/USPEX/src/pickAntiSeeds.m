function pickAntiSeeds()
global ANTISEEDS
global ORG_STRUC
ANTISEEDS = struct('FINGERPRINT',{},'Sigma',{},'Max',{});
ANTISEEDS(1).FINGERPRINT = [];
ANTISEEDS(1).Sigma = [];
ANTISEEDS(1).Max = [];
cd AntiSeeds
disp('Read AntiSeeds ...')
try
[fid,message] = fopen('POSCARS');
loops_seed = 0;
while 1
loops_seed = loops_seed + 1;
tmp = fgetl(fid); 
if tmp == -1; break; end
scale_factor = fgetl(fid);           
optlattice = fscanf(fid,'%g',[3,3]); % lattice vectors
optlattice = optlattice'*str2num(scale_factor);
tmp = fgetl(fid);                    
atomType = fgetl(fid);   
ntyp = fgetl(fid);       
ntyp = str2num(ntyp);
natom = sum(ntyp);       
tmp = fgetl(fid);        
candidate_pop = fscanf(fid,'%g',[3,natom]);
candidate_pop = candidate_pop';
tmp = fgetl(fid);
numIons = ntyp;
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(optlattice, candidate_pop, numIons, GetElement(length(ntyp), atomType));
[order, fing, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
ANTISEEDS(loops_seed).FINGERPRINT = fing;
ANTISEEDS(loops_seed).Sigma = ORG_STRUC.antiSeedsSigma;
ANTISEEDS(loops_seed).Max   = ORG_STRUC.antiSeedsMax;
disp([' --> AntiSeed ' num2str(loops_seed) ' added ...'])
end
status = fclose(fid);
catch
end
fclose('all');
cd (ORG_STRUC.homePath)
safesave ('ANTISEEDS.mat', ANTISEEDS);
