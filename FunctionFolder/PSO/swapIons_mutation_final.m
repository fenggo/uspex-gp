function [goodSon, swapNow] = swapIons_mutation_final(father,numIons,order)
global ORG_STRUC
ordering_on = 0; 
ionChange = zeros(1,length(numIons));
for ind = 1:length(numIons)
ionChange(ind) = sum(numIons(1:ind));
end
if ORG_STRUC.howManySwaps < 1
ORG_STRUC.howManySwaps = 1;
end
swapNow = RandInt(1,1,[1,ORG_STRUC.howManySwaps]);
for any = 1:swapNow
if ORG_STRUC.specificSwaps==0
firstIon = RandInt(1,1,[1,sum(numIons)]);
whichTypeFirst = find (ionChange+0.5>firstIon);
else
specified = 1;
while specified
firstIon = RandInt(1,1,[1,sum(numIons)]);
whichTypeFirst = find (ionChange+0.5>firstIon);
if(~isempty(whichTypeFirst))
if ~isempty(find(whichTypeFirst(1) == ORG_STRUC.specificSwaps))
[row,col]=find(whichTypeFirst(1) == ORG_STRUC.specificSwaps);
acceptedSwaps=[];
for counting =1:length(row)
if  col(counting)-0.1>1
acceptedSwaps=cat(1,acceptedSwaps,ORG_STRUC.specificSwaps(row(counting),1:col(counting)-1));
end
if col(counting)+0.1<size(ORG_STRUC.specificSwaps,2)
acceptedSwaps=cat(1,acceptedSwaps,ORG_STRUC.specificSwaps(row(counting),col(counting)+1:end));
end
end
specified=0;
end
end
end
end
match = 1;
while match
secondIon = RandInt(1,1,[1,sum(numIons)]);
whichTypeSecond = find (ionChange+0.5>secondIon);
if numIons(whichTypeFirst(1)) == sum(numIons)
match = 0;
end;
if ORG_STRUC.specificSwaps == 0
if whichTypeFirst(1) ~= whichTypeSecond(1)
match = 0;
end
else
if ~isempty(find(whichTypeSecond(1)==acceptedSwaps))
match = 0;
end
end
end
if ordering_on
tmp = zeros(numIons(whichTypeFirst(1)),1);
for i = 1 : numIons(whichTypeFirst(1))
tmp(i) = sum(numIons(1:whichTypeFirst(1)))-i+1;
end
atoms = sort_order(order, tmp, 'ascend');
firstIon = tmp(atoms(1));
tmp = zeros(numIons(whichTypeSecond(1)),1);
for i = 1 : numIons(whichTypeSecond(1))
tmp(i) = sum(numIons(1:whichTypeSecond(1)))-i+1;
end
atoms = sort_order(order, tmp, 'ascend');
secondIon = tmp(atoms(1));
end
tempFather = father(firstIon,:) ;
father(firstIon,:) = father(secondIon,:);
father(secondIon,:) = tempFather;
end
goodSon = father;
