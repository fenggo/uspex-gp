function [reduced, newCoord, newLat, newNumIons, newSupercell] = reduceSupercell(lat, coord, numIons, supercell, tolerance)
maxMultiplicity = 20; 
reduced = 0; 
N = sum(numIons);
atomType = zeros(1,N);
counter = 1;
for i = 1 : size(numIons)
for j = 1 : numIons(i)
atomType(counter) = i;
counter = counter + 1;
end
end
try
Mult = [1 1 1]; 
for Axis = 1 : 3 
for m = 2 : maxMultiplicity 
if mod(supercell(Axis), m) > 0
continue
end
pairFound = zeros(1,N);
for a1 = 1 : N
vec_a1 = coord(a1,:);
vec_a1(Axis) = vec_a1(Axis) + (1/m);
vec_a1 = vec_a1 - floor(vec_a1); 
for a2 = 1 : N
diff_vec = (coord(a2,:) - vec_a1)*lat;       
if (norm(diff_vec) <= tolerance) & (atomType(a1) == atomType(a2)) 
pairFound(a1) = 1;
break;
end
end
end
if (sum(pairFound) == N)
Mult(Axis) = m;  
end
end
end
if Mult(1)*Mult(2)*Mult(3) > 1
reduced = 1;
newLat(1,:) = lat(1,:)/Mult(1);
newLat(2,:) = lat(2,:)/Mult(2);
newLat(3,:) = lat(3,:)/Mult(3);
newSupercell(1) = supercell(1)/Mult(1);
newSupercell(2) = supercell(2)/Mult(2);
newSupercell(3) = supercell(3)/Mult(3);
newNumIons = numIons/(Mult(1)*Mult(2)*Mult(3));
newCoord = zeros(sum(newNumIons),3); 
counter = 1;
c = coord - floor(coord); 
for a = 1 : N
if (c(a,1) < 1/Mult(1)) & (c(a,2) < 1/Mult(2)) & (c(a,3) < 1/Mult(3))
newCoord(counter,:) = c(a,:);
newCoord(counter,1) = newCoord(counter,1) * Mult(1);
newCoord(counter,2) = newCoord(counter,2) * Mult(2);
newCoord(counter,3) = newCoord(counter,3) * Mult(3);
counter = counter + 1;
end
end
if sum(newNumIons) ~= size(newCoord,1)
reduced = 0; 
end
else
newCoord = coord;
newLat = lat;
newNumIons = numIons;
newSupercell = supercell;
end
catch
reduced = 0;
newCoord = coord;
newLat = lat;
newNumIons = numIons;
newSupercell = supercell;
end