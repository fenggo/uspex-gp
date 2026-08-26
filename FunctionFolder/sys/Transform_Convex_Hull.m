function convex_hull0 = Transform_Convex_Hull(convex_hull, numIons, flag);
N_Block = size(numIons, 1);
N_Atom  = size(numIons, 2);
N_Point = size(convex_hull, 1);
if flag == 0
convex_hull0 = zeros(N_Point, N_Atom+2);
for i = 1 : N_Point
Block_tmp = convex_hull(i, 1:N_Block);
Atom_tmp  = Block_tmp*numIons;
convex_hull0(i,1:N_Atom) = Atom_tmp;
convex_hull0(i,N_Atom+1) = convex_hull(i, N_Block+1)*sum(Block_tmp)/sum(Atom_tmp);  
convex_hull0(i,N_Atom+2) = convex_hull(i, N_Block+2);
end
else 
convex_hull0 = zeros(N_Point, N_Block+2);
for i = 1 : N_Point
Atom_tmp = convex_hull(i, 1:N_Atom);
Block_tmp = Atom_tmp/numIons;
convex_hull0(i,1:N_Block) = Block_tmp;
convex_hull0(i,N_Block+1) = convex_hull(i, N_Atom+1)*sum(Atom_tmp)/sum(Block_tmp);  
convex_hull0(i,N_Block+2) = convex_hull(i, N_Atom+2);
end
end
