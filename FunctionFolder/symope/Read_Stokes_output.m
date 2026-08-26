function [candidate, lattice, errorS] = Read_Stokes_output(filename, numIons)
candidate = zeros(sum(numIons),3);
lattice = [1 1 1 0 0 0]';
errorS = 0;
if exist(filename)
try
handle = fopen(filename);
f1 = ceil(length(numIons)^2/10);
for LL = 1 : 6 + f1
tmp = fgetl(handle); 
end
lat_tmp = fgetl(handle); 
error1 = findstr(lat_tmp, 'error');
if isempty(error1)
lattice = sscanf(lat_tmp,'%*s %*s %g %g %g %g %g %g');
tmp = fgetl(handle); 
tmp = fgetl(handle); 
for LL = 1 : sum(numIons)
atom_tmp = fgetl(handle); 
atom_1 = sscanf(atom_tmp,'%*s %*s %g %g %g');
candidate(LL,1) = atom_1(1); candidate(LL,2) = atom_1(2); candidate(LL,3) = atom_1(3);
end
fclose(handle);
else
disp(['Error while generating crystal with symmetry ' num2str(nsym) ' and N_atoms ' num2str(numIons)]);
fclose(handle);
errorS = 1;
end
catch
fclose(handle);
errorS = 1;
end
else
disp(['ooopps, ' filename 'does not exist']);
errorS = 1;
end
