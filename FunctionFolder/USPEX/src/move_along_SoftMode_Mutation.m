function [lat, new_Coord, deviation]= move_along_SoftMode_Mutation(coord, numIons, lattice, eigvector, mut_degree)
global ORG_STRUC
global POP_STRUC
N = sum(numIons);
vec = zeros(1,3);
new_Coord = coord;
lat = lattice;
close_enough = -0.37*log(ORG_STRUC.goodBonds);        
vect = zeros(1,3);
R_val = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
s = covalentRadius(ceil(ORG_STRUC.atomType(i)));
R_val(i) = str2num(s);
end
if ORG_STRUC.dimension==0 
coord = coord + 0.0001;
[lat, candidate] = makeCluster(lattice, coord, ORG_STRUC.vacuumSize(1));
[candidates] = moveCluster(lat, candidate);
lattice = lat;
coord = candidates;
end
at_types = zeros(1,N);
for k = 1 : N
tmp = k;
while tmp > 0
at_types(k) = at_types(k) + 1;
tmp = tmp - numIons(at_types(k));
end
end                                                  
if length(lattice) == 6
lattice = latConverter(lattice);
end
coef = 0;
for i = 1 : N
vec(1) = eigvector((i-1)*3+1);
vec(2) = eigvector((i-1)*3+2);
vec(3) = eigvector((i-1)*3+3);
if norm(vec) > coef
coef = norm(vec);
end
end
normfac = mut_degree*ORG_STRUC.howManyMut/coef;   
notDone = 1;
step = 0;
while notDone
for i = 1 : N                      
vec(1) = eigvector((i-1)*3+1);
vec(2) = eigvector((i-1)*3+2);
vec(3) = eigvector((i-1)*3+3);
new_Coord(i,:) = coord(i,:) + normfac*vec*inv(lattice);
deviation(i) = norm(normfac*vec);
end
if ORG_STRUC.dimension==0 
badAtoms = ones(1,N);
dist = 10*ones(N,N);
min_dist = zeros(1,N);
for i = 1 : N
for j = i+1 : N 
vect(1) = new_Coord(i,1) - new_Coord(j,1);
vect(2) = new_Coord(i,2) - new_Coord(j,2);
vect(3) = new_Coord(i,3) - new_Coord(j,3);
delta = sqrt(sum((vect*lattice).^2)) - R_val(at_types(i)) - R_val(at_types(j));
dist(i,j) = delta;
dist(j,i) = delta;
if delta < close_enough(at_types(i), at_types(j))
badAtoms(i) = 0;
badAtoms(j) = 0;
end 
end
min_dist(i) = min(dist(i,:));
end
if sum(badAtoms) == 0
notDone = 0;
else
step = step + 1;
for i = 1 : N
if badAtoms(i)
eigvector((i-1)*3+1) = eigvector((i-1)*3+1)*0.9;
eigvector((i-1)*3+2) = eigvector((i-1)*3+2)*0.9;
eigvector((i-1)*3+3) = eigvector((i-1)*3+3)*0.9;
end
end
end
else
notDone = 0;
end
if step > 10
%unix(['echo main_structure >> POSCAR_cluster']);
%unix(['echo ' '1.0' ' >> POSCAR_cluster']);
%  unix(['echo ' num2str(lattice(latticeLoop,:)) ' >> POSCAR_cluster']);
%unix(['echo ' num2str(numIons) ' >> POSCAR_cluster']);
%unix(['echo ' 'Direct' ' >> POSCAR_cluster']);
%  unix(['echo ' num2str(coord(coordLoop,:)) ' >> POSCAR_cluster']);
notDone = 0;
end
end
