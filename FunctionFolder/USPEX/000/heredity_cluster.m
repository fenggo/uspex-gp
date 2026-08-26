function [numIons, offspring,potentialLattice,fracFrac,dimension,offset,fracLattice] = heredity_cluster(parent1,parent2,parentlat1,parentlat2, order1, order2, fracFrac)
global POP_STRUC
global ORG_STRUC
ordering_on = ORG_STRUC.ordering;
numIons = ORG_STRUC.numIons;
fracFrac = 0.49999999;
fracLattice = fracFrac;
dimension = RandInt(1,1,[1,3]);
offset = pi*(rand(3,1) - 0.5);
[parentlat1, parent1] = rotateCluster(parent1, parentlat1, 0, 0, offset(1,1));
[parentlat2, parent2] = rotateCluster(parent2, parentlat2, 0, 0, offset(1,1));
if rand(1) > ORG_STRUC.percSliceShift   
else 
[parentlat1, parent1] = rotateCluster(parent1, parentlat1, offset(2,1), offset(3,1), 0);
[parentlat2, parent2] = rotateCluster(parent2, parentlat2, offset(2,1), offset(3,1), 0);
end
whichIons1 = find(parent1(:,dimension)<fracFrac);
whichIons2 = find(parent2(:,dimension)>=fracFrac);
correlation_coefficient = ORG_STRUC.correlation_coefficient;
cor_dir = ORG_STRUC.cor_dir;
ord_parent = rand; 
if ord_parent < 0.5 
nAt = size(parent1,1);
Lchar = 0.5*power(abs(det(parentlat1))/nAt, 1/3); 
else
nAt = size(parent2,1);
Lchar = 0.5*power(abs(det(parentlat2))/nAt, 1/3); 
end
Nmax = 0.5*round((6*pi*pi*nAt)^1/3);
Nslabs = round(Nmax/(1+(Nmax-1)*(cos(correlation_coefficient*pi/2))^2));
%disp(['Number of cuts = ' num2str(Nslabs)])
if ordering_on
if length(whichIons1) > 0
ord1 = sum(order1(whichIons1))/length(whichIons1);
elseif cor_dir <= 0
ord1 = 0;
else
ord1 = 1;
end
if length(whichIons2) > 0
ord2 = sum(order2(whichIons2))/length(whichIons2);
elseif cor_dir <= 0
ord2 = 0;
else
ord2 = 1;
end
lat1_candidate = parentlat1;
lat2_candidate = parentlat2;
parent1_candidate = parent1;
parent2_candidate = parent2;
for i = 1 : Nslabs
offset_extra = pi*(rand(3,1) - 0.5);
[lat1_tmp, parent1_tmp] = rotateCluster(parent1, parentlat1, offset_extra(1,1),offset_extra(2,1),offset_extra(3,1));
[lat2_tmp, parent2_tmp] = rotateCluster(parent2, parentlat2, offset_extra(1,1),offset_extra(2,1),offset_extra(3,1));
whichIons1_extra = find(parent1_tmp(:,dimension)<fracFrac);
whichIons2_extra = find(parent2_tmp(:,dimension)>=fracFrac);
if length(whichIons1_extra) > 0
ord1_extra = sum(order1(whichIons1_extra))/length(whichIons1_extra);
elseif cor_dir <= 0
ord1_extra = 0;
else
ord1_extra = 1;
end
if length(whichIons2_extra) > 0
ord2_extra = sum(order2(whichIons2_extra))/length(whichIons2_extra);
elseif cor_dir <= 0
ord2_extra = 0;
else
ord2_extra = 1;
end
if ord_parent < 0.5 
if ((ord1 > ord1_extra) & (cor_dir <= 0)) | ((ord1 < ord1_extra) & (cor_dir == 1))
else
whichIons1 = whichIons1_extra;
whichIons2 = whichIons2_extra;
ord1 = ord1_extra;
ord2 = ord2_extra;
lat1_candidate = lat1_tmp;
lat2_candidate = lat2_tmp;
parent1_candidate = parent1_tmp;
parent2_candidate = parent2_tmp;
end
else 
if ((ord2 > ord2_extra) & (cor_dir <= 0)) | ((ord2 < ord2_extra) & (cor_dir == 1))
else
whichIons1 = whichIons1_extra;
whichIons2 = whichIons2_extra;
ord1 = ord1_extra;
ord2 = ord2_extra;
lat1_candidate = lat1_tmp;
lat2_candidate = lat2_tmp;
parent1_candidate = parent1_tmp;
parent2_candidate = parent2_tmp;
end
end
end
parentlat1 = lat1_candidate;
parentlat2 = lat2_candidate;
parent1 = parent1_candidate;
parent2 = parent2_candidate;
end
offspring = zeros(0,3);
ionCount = zeros(length(ORG_STRUC.atomType),1);
for ind = 1 : length(numIons)
ionChange(ind) = sum(numIons(1:ind));
end
ionCh1 = zeros(1,length(ionChange)+1);
ionCh1(2:end) = ionChange;
ionCh2 = ionCh1;
for ind = 1:length(ORG_STRUC.atomType)
candidates = [];
smaller1 = find(whichIons1<=ionCh1(ind+1));
smaller2 = find(whichIons2<=ionCh2(ind+1));
ions1  = find(whichIons1(smaller1)>ionCh1(ind));
ions2  = find(whichIons2(smaller2)>ionCh2(ind));
ionCount(ind,1) = length(ions1) + length(ions2);
if ionCount(ind,1) < numIons(ind)
howmany = numIons(ind) - ionCount(ind,1);
suppleIons_1 = length(find(rand(howmany,1)<fracFrac));
suppleIons_2 = howmany-suppleIons_1;
if suppleIons_1 > 0
tempInd = find(parent2(ionCh1(ind)+1:ionCh1(ind+1),dimension)<fracFrac);
if cor_dir <= 0  
atoms = sort_order(order2, tempInd, 'descend');
else 
atoms = sort_order(order2, tempInd, 'ascend');
end
candidates = parent2(ionCh2(ind) + tempInd(atoms(1:suppleIons_1)) , :);
end
if suppleIons_2 > 0
tempInd = find(parent1(ionCh1(ind)+1:ionCh1(ind+1),dimension)>fracFrac);
if cor_dir <= 0  
atoms = sort_order(order1, tempInd, 'descend');
else 
atoms = sort_order(order1, tempInd, 'ascend');
end
candidates = cat(1,candidates,parent1(ionCh1(ind) + tempInd(atoms(1:suppleIons_2)) , :));
end
elseif  ionCount(ind,1)>numIons(ind)
howmany = ionCount(ind,1) - numIons(ind);
deleteIons_1 = length(find(rand(howmany,1)<fracFrac));
deleteIons_2 = howmany - deleteIons_1;
if length(ions1) < deleteIons_1
oops = deleteIons_1 - length(ions1);
deleteIons_1 = length(ions1);
deleteIons_2 = deleteIons_2+oops;
elseif length(ions2) < deleteIons_2
oops = deleteIons_2 - length(ions2);
deleteIons_2 = length(ions2);
deleteIons_1 = deleteIons_1+oops;
end
for xy = 1:deleteIons_1
ions1(RandInt(1,1,[1,length(ions1)])) = [];
end
for xy = 1:deleteIons_2
ions2(RandInt(1,1,[1,length(ions2)])) = [];
end
end
addOn = cat(1,parent1(whichIons1(smaller1(ions1(:))) , :), parent2(whichIons2(smaller2(ions2(:))) , :));
addOn = cat(1,addOn,candidates);
offspring = cat(1,offspring,addOn);
end
temp_potLat = fracFrac*parentlat1 + (1-fracFrac)*parentlat2;
volLat = det(temp_potLat);
latVol = det(parentlat1)*fracFrac + det(parentlat2)*(1-fracFrac);
ratio = latVol/volLat;
temp_potLat = latConverter(temp_potLat);
temp_potLat(1:3)= temp_potLat(1:3)*(ratio)^(1/3);
potentialLattice = latConverter(temp_potLat);
