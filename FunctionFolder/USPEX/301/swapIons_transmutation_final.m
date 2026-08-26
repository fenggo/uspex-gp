function [goodSon, goodFather, transNow, numIons, numBlocks] = swapIons_transmutation_final(Ind_No)
global POOL_STRUC
global ORG_STRUC
numBlocks = POOL_STRUC.POPULATION(Ind_No).numBlocks;
addBlocks = zeros(1, length(numBlocks));  
delBlocks = addBlocks;                    
father = POOL_STRUC.POPULATION(Ind_No).COORDINATES;
goodFather = father;
numIons = POOL_STRUC.POPULATION(Ind_No).numIons;
howManyTrans = round(ORG_STRUC.howManyTrans*sum(numBlocks));
if howManyTrans < 1 | isnan( howManyTrans )
howManyTrans = 1;
end
transNow = RandInt(1,1,[1,howManyTrans]);
if ~(ORG_STRUC.specificTrans == 0)
validStr = 0;
for i = 1 : length(numBlocks)
if (~isempty(find(i == ORG_STRUC.specificTrans))) & (numBlocks(i) > 0)
validStr = 1;
end
end
if validStr == 0
transNow = 0;
goodSon = [];
end
end
for any = 1 : transNow
if ORG_STRUC.specificTrans == 0
constrain = 1;
while constrain  
whichTypeFirst = RandInt(1,1,[1,length(numBlocks)]);
if numBlocks(whichTypeFirst) > 0
constrain = 0;
end
end
else
specified = 1;
while specified
whichTypeFirst = RandInt(1,1,[1,length(numBlocks)]);
if (~isempty(find(whichTypeFirst == ORG_STRUC.specificTrans))) & (numBlocks(whichTypeFirst) > 0)
[row, col] = find(whichTypeFirst == ORG_STRUC.specificTrans);
acceptedTrans = [];
for counting = 1 : length(row)
if col(counting)-0.1 > 1
acceptedTrans = cat(1,acceptedTrans,ORG_STRUC.specificTrans(row(counting),1:col(counting)-1));
end
if col(counting)+0.1 < size(ORG_STRUC.specificTrans,2)
acceptedTrans = cat(1,acceptedTrans,ORG_STRUC.specificTrans(row(counting),col(counting)+1:end));
end
end
specified = 0;
end
end
end
match = 1;
while match
whichTypeSecond = RandInt(1,1,[1,length(numBlocks)]);
if ORG_STRUC.specificTrans == 0
if (whichTypeFirst ~= whichTypeSecond) 
match = 0;
end
else
if ~isempty(find(whichTypeSecond == acceptedTrans))
match = 0;
end
end
end
delBlocks(whichTypeFirst) = delBlocks(whichTypeFirst) + 1;
addBlocks(whichTypeSecond) = addBlocks(whichTypeSecond) + 1;
numBlocks(whichTypeFirst) = numBlocks(whichTypeFirst) - 1;
numBlocks(whichTypeSecond) = numBlocks(whichTypeSecond) + 1;
end
if transNow > 0   
addAtoms = addBlocks*ORG_STRUC.numIons; 
delAtoms = delBlocks*ORG_STRUC.numIons; 
for i = 1 : length(numIons)
if delAtoms(i) > numIons(i)
addAtoms(i) = addAtoms(i) - (delAtoms(i) - numIons(i));
delAtoms(i) = numIons(i);
end
end
toDelete = zeros(1, sum(delAtoms));  
toAdd = zeros(1, sum(addAtoms));     
atomTypes = zeros(1, sum(numIons));  
c = 0;
for i = 1 : length(numIons)
dummy = randperm(numIons(i)) + c;
c = c + numIons(i);
i1 = sum(delAtoms(1:i-1));
toDelete(i1+1 : i1+delAtoms(i)) = dummy(1:delAtoms(i));
i1 = sum(addAtoms(1:i-1));
toAdd(i1+1 : i1+addAtoms(i)) = i;
i1 = sum(numIons(1:i-1));
atomTypes(i1+1 : i1+numIons(i)) = i;
end
atomTypes(toDelete(1:end)) = 0;
dummy = randperm(length(toAdd));
i1 = 1;
for i = 1 : sum(numIons)
if (atomTypes(i) == 0) & (i1 <= length(toAdd))
atomTypes(i) = toAdd(dummy(i1));
i1 = i1 + 1;
end
end
numIons = numBlocks*ORG_STRUC.numIons;   
goodSon = zeros(sum(numIons), 3);
g = 1;
if length(toDelete) >= length(toAdd)    
for i = 1 : length(numIons)
for j = 1 : length(atomTypes)
if atomTypes(j) == i
goodSon(g, 1:3) = father(j, 1:3);
g = g + 1;
end
end
end
else   
extraAtoms = rand(length(toAdd) - length(toDelete),3);
extraTypes = toAdd(dummy(i1:end)); 
for i = 1 : length(numIons)
for j = 1 : length(atomTypes)
if atomTypes(j) == i
goodSon(g, 1:3) = father(j, 1:3);
g = g + 1;
end
end
for j = 1 : length(extraTypes)
if extraTypes(j) == i
goodSon(g, 1:3) = extraAtoms(j, 1:3);
g = g + 1;
end
end
end
end
end
