function [coor, lat] = Read_QE_Structure()
bohr = 0.529177; 
[nothing, ans] = unix('grep CELL_P output |tail -1');
if ~isempty(findstr('alat',ans))  
tmp = ans(23:32);
scale_fac = str2num(tmp);
lat = callAWK('QE_cell.awk', '| tail -3', 'output');
else  
[nothing, scaleStr] = unix('./getStuff output "lattice parameter" 6');
scale_fac = str2num(scaleStr);
lat = callAWK('QE_cell.awk', '| tail -3', 'output');
end
lat = lat*bohr;
lat = lat*scale_fac;
[a,SCF]=unix('grep ATOMIC_POSITIONS output |wc -l');  
scf = str2num(SCF);
COOR = callAWK('QE_atom.awk','output');
numIons = size(COOR,1)/scf;
if scf > 1
COOR(1:numIons*(scf-1),:)=[];
end
coor = COOR;
