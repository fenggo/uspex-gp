function [resulted_sg_no, resulted_sg_name] = ...
determine_spacegroup(lattice, coordinates, atom_type, ...
num_atoms, nsym, structure_id, ...
tolerance, home_path, method)
if nargin < 9
method = 'Stokes';  
end
error_phonopy = 1;
if isequal(method, 'Phonopy')
[error_phonopy, nothing] = unix('which phonopy');
end
if error_phonopy == 0
structures_dir = [home_path '/structures'];
if ~isequal(exist(structures_dir, 'dir'), 7)  
[nothing, nothing] = unix(['mkdir ' structures_dir]);
end
structure_id = num2str(structure_id, '%04i');
outfile = [structures_dir '/' structure_id '_SG_' num2str(nsym)];  % output POSCAR file
atoms = '';
for i=1:size(atom_type, 2)
if i == 1
atoms = [char(megaDoof(atom_type(i)))];
else
atoms = [atoms ' ' char(megaDoof(atom_type(i)))];
end
end
atoms = {atoms};
content = generate_poscar(lattice, coordinates, atoms, num_atoms, outfile, nsym);
phonopy_command = ['ulimit -s unlimited; phonopy --symmetry --tolerance=' num2str(tolerance) ...
' -c ' outfile ' | grep space_group_type | cut -d: -f2- | awk ''{print $1,$2}'' '];
[nothing, resulted_sg] = unix(phonopy_command);
try
parsed = sscanf(resulted_sg, '%s (%i)');
resulted_sg_no   = parsed(end);
resulted_sg_name = char(parsed(1:end-1)');
catch
resulted_sg_no   = 1;
resulted_sg_name = 'P1';
end
else
method = 'Stokes';
resulted_sg_no = findsym_stokes(lattice, coordinates, num_atoms, atom_type, tolerance);
if ~isempty(resulted_sg_no)
resulted_sg_no = str2num(resulted_sg_no);
else
resulted_sg_no = 1;
end
resulted_sg_name = spaceGroups(resulted_sg_no);
end
% disp(' ');
disp(['  Initial space group: ' spaceGroups(nsym) ' (' num2str(nsym)           ')']);
disp(['  Actual  space group: ' resulted_sg_name  ' (' num2str(resulted_sg_no) ...
') (determined with tolerance=' num2str(tolerance) ')']);
