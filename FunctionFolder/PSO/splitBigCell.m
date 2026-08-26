function [lat, errorS, coordinates] = splitBigCell(startLat, splitInto, numIons, nsym)
errorS = 0;
global ORG_STRUC
global POP_STRUC
splitByCellNumber = 1;
if ORG_STRUC.constLattice
startLat = latConverter(ORG_STRUC.lattice);
end
b = numIons;
L = length(b);
b = sort(b);
tmp = b;
for i = 1:length(b)
b(end-i+1) = tmp(i);
end
while sum(b(2:L)) ~= 0
[m,k] = max(b(2:L));
if b(1) > m
b(1) = b(1) - m;
else
b(k+1) = b(k+1) - b(1);
end
end
GrCD = b(1);
b = numIons;
nStrEnt = GrCD;
nFull = sum(numIons); 
nAt = nFull/GrCD; 
c = b/GrCD;
[lattice, IX] = sort(startLat(1:3)); 
opt = 10*nFull;
opt1 = 10*nFull;
x = 0;
y = 0;
z = 0;
if splitByCellNumber == 1
i1 = splitInto.^(1/3);
i2 = splitInto.^(0.5);
coor = zeros(3,1);
for i = 1 : i1+1
for j = i : i2+1
k = floor(splitInto/(i*j));
if k >= j
candidate = abs(splitInto-i*j*k);
split_lattice_match = i/lattice(1)+lattice(1)/i+j/lattice(2)+lattice(2)/j+k/lattice(3)+lattice(3)/k;
if (candidate < opt) | ((candidate == opt) & (split_lattice_match < opt1))
coor(1) = i; coor(2) = j; coor(3) = k; opt = candidate; opt1 = split_lattice_match;
end
end
end
end
[l1, IX1]=sort(IX);
x = coor(IX1(1));
y = coor(IX1(2));
z = coor(IX1(3));
else   
end
doAtoms = zeros(1,L); 
addAtoms = zeros(1,L); 
genAtoms = zeros(1,L); 
for atomType = 1 : L
doAtoms(atomType) = floor(numIons(atomType)/splitInto);
addAtoms(atomType) = numIons(atomType) - splitInto*doAtoms(atomType);
if addAtoms(atomType) > 0
genAtoms(atomType) = doAtoms(atomType) + 1; 
else
genAtoms(atomType) = doAtoms(atomType);
end
end
coord_splitter = rand(sum(genAtoms(:)),3);
if ORG_STRUC.nsymN(1,1) == 0  
lattice1 = startLat;
lattice1(1) = lattice1(1)/x;
lattice1(2) = lattice1(2)/y;
lattice1(3) = lattice1(3)/z;
for l1 = 1 : 5
cd([ORG_STRUC.homePath '/CalcFoldTemp'])
[coord_splitter, lat1, errorS] = symope_crystal(nsym, genAtoms, lattice1, ORG_STRUC.minDistMatrice, ORG_STRUC.sym_coef);
cd(ORG_STRUC.homePath)
if errorS == 0
break
end
end
if errorS == 1
coord_splitter = rand(sum(genAtoms(:)),3);
elseif ORG_STRUC.constLattice == 0
lat1 = latConverter(lat1);
startLat(1) = lat1(1)*x;  
startLat(2) = lat1(2)*y;
startLat(3) = lat1(3)*z;
startLat(4:6) = lat1(4:6);
end
end
coordinates = rand(nFull,3);
counter = 1;
counter1 = 0;
for atomType = 1 : L
for a1 = 1 : doAtoms(atomType)
counter1 = counter1 + 1;
for i = 1 : x
for j = 1 : y
for k = 1 : z
coordinates(counter,1) = (i-1+coord_splitter(counter1,1))/x;
coordinates(counter,2) = (j-1+coord_splitter(counter1,2))/y;
coordinates(counter,3) = (k-1+coord_splitter(counter1,3))/z;
counter = counter + 1;
end
end
end
end
p = randperm(splitInto);
if addAtoms(atomType) > 0
counter1 = counter1 + 1;
end
for a1 = 1 : addAtoms(atomType)
div1 = floor(p(a1)/(x*y));
mod1 = p(a1) - div1*x*y;
if mod1 == 0
k = div1; p(a1) = x*y;
else
k = div1 + 1; p(a1) = mod1;
end;
div1 = floor(p(a1)/x);
mod1 = p(a1) - div1*x;
if mod1 == 0
j = div1; i = x;
else
j = div1 + 1; i = mod1;
end;
coordinates(counter,1) = (i-1+coord_splitter(counter1,1))/x;
coordinates(counter,2) = (j-1+coord_splitter(counter1,2))/y;
coordinates(counter,3) = (k-1+coord_splitter(counter1,3))/z;
counter = counter + 1;
end
end
lat = latConverter(startLat);
if errorS == 0
disp(['split into: x = ' num2str(x) ', y = ' num2str(y) ', z = ' num2str(z)])
end
