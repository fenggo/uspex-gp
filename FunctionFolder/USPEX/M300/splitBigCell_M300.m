function [coordinates, lat, split, errorS] = splitBigCell_M300(startLat, splitInto, numIons, nsym, minDistMatrice)
lat = startLat;
lat_6 = latConverter(startLat);
N_type = length(numIons);
GrCD = numIons(1);
if N_type > 1
for i = 2:N_type
GrCD = gcd(GrCD, numIons(i));
end
end
nFull = sum(numIons); 
nAt = nFull/GrCD; 
[lattice, IX] = sort(lat_6(1:3)); 
opt = 10*nFull;
opt1 = 10*nFull;
x=0; y=0; z=0;
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
split(1) = coor(IX1(1));  split(2)=coor(IX1(2)); split(3)=coor(IX1(3));
doAtoms =  zeros(1,N_type); 
addAtoms = zeros(1,N_type); 
genAtoms = zeros(1,N_type); 
for atomType = 1 : N_type
doAtoms(atomType) = floor(numIons(atomType)/splitInto);
addAtoms(atomType) = numIons(atomType) - splitInto*doAtoms(atomType);
if addAtoms(atomType) > 0
genAtoms(atomType) = doAtoms(atomType) + 1; 
else
genAtoms(atomType) = doAtoms(atomType);
end
end
lattice1 = startLat;
lattice1(1) = lattice1(1)/split(1); lattice1(2) = lattice1(2)/split(2); lattice1(3) = lattice1(3)/split(3);
cd(['CalcFoldTemp'])
[coord_splitter, lat1, errorS] = symope_GB(nsym, genAtoms, lattice1, minDistMatrice);
cd ..
coordinates = rand(nFull,3);
if errorS ==0
counter = 1;
counter1 = 0;
for atomType = 1 : N_type
for a1 = 1 : doAtoms(atomType)
counter1 = counter1 + 1;
for i = 1 : split(1)
for j = 1 : split(2)
for k = 1 : split(3)
coordinates(counter,1) = (i-1+coord_splitter(counter1,1))/split(1);
coordinates(counter,2) = (j-1+coord_splitter(counter1,2))/split(2);
coordinates(counter,3) = (k-1+coord_splitter(counter1,3))/split(3);
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
div1 = floor(p(a1)/(split(1)*split(2)));
mod1 = p(a1) - div1*split(1)*split(2);
if mod1 == 0
k = div1; p(a1) = split(1)*split(2);
else
k = div1 + 1; p(a1) = mod1;
end;
div1 = floor(p(a1)/split(1));
mod1 = p(a1) - div1*split(1);
if mod1 == 0
j = div1; i = split(1);
else
j = div1 + 1; i = mod1;
end;
coordinates(counter,1) = (i-1+coord_splitter(counter1,1))/split(1);
coordinates(counter,2) = (j-1+coord_splitter(counter1,2))/split(2);
coordinates(counter,3) = (k-1+coord_splitter(counter1,3))/split(3);
counter = counter + 1;
end
end
end
