function [lat_conv, coor_conv] = unitCellFromPrimitive(cellType, lat_prim, coor_prim)
if cellType == 'F' 
U = [-1 +1 +1; +1 -1 +1; +1 +1 -1];
Duplicate = [0 0 0; 0.5 0 0.5; 0 0.5 0.5; 0.5 0.5 0];
elseif cellType == 'I' 
U = [0 1 1; 1 0 1; 1 1 0];
Duplicate = [0 0 0; 0.5 0.5 0.5];
elseif cellType == 'C'
U = [1 1 0; 1 -1 0; 0 0 1];
Duplicate = [0 0 0; 0.5 0.5 0];
elseif cellType == 'B' 
U = [0 0 1; 1 -1 0; 1 1 0];
Duplicate = [0 0 0; 0 0.5 0.5];
elseif cellType == 'A'
U = [0 1 -1; 0 1 1; 1 0 0];
Duplicate = [0 0 0; 0.5 0.5 0];
elseif cellType == 'R' 
U = [1 -1 0; 0 1 -1; 1 1 1];
Duplicate = [0 0 0; -1/3 1/3 1/3; -2/3 2/3 2/3];
end
N_atom = size(coor_prim,1);
N_Dupl = size(Duplicate, 1);
lat_conv = U*lat_prim;    
coor_prim = coor_prim*lat_prim/lat_conv;  
coor_conv = zeros(N_atom*N_Dupl, 3);
for i = 1: N_atom
for j = 1: N_Dupl
coor_conv(N_Dupl*(i-1)+j,:) = coor_prim(i,:) + Duplicate(j,:);
end
end
coor_conv = coor_conv - floor(coor_conv);
lat_conv = latConverter(latConverter(lat_conv)); 
