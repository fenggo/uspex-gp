function [numIons, numBlocks] = findDesiredComposition(maxBlocks, blocks, child)
maxAtoms = maxBlocks*blocks;
maxAdded = maxAtoms - child; 
Nb = size(blocks,1); 
Nt = size(blocks,2); 
tmp = 1; 
for i = 1 : Nb
min1 = max(max(maxAtoms));
for j = 1 : Nt
if blocks(i,j) > 0
if min1 > maxAtoms(j)/blocks(i,j)
min1 = floor(maxAtoms(j)/blocks(i,j));
end
end
end
if min1 > 0
tmp = tmp*min1;
end
end
if tmp > 0   
bestGreed = sum(maxAtoms);
for g = 1 : 20
tolerance = round(rand*2);    
blockN = zeros(1,Nb); 
child_tmp = child;
crutch = rand(Nb,1);
[junk, blockOrder] = sort(crutch);
for i = 1 : Nb
block = blocks(blockOrder(i), 1:end);
ind = 1;
min1 = sum(child_tmp);
while 1
if min(child_tmp - ind*block) < -1*tolerance
break;
end
if min(maxAdded + (child_tmp - ind*block)) < 0
break;
end
if min1 > sum(abs(child_tmp - ind*block))
min1 = sum(abs(child_tmp - ind*block));
blockN(blockOrder(i)) = ind;
end
ind = ind + 1;
end
child_tmp = child_tmp - blockN(blockOrder(i))*block;
end
if bestGreed > sum(abs(child - blockN*blocks))
bestGreed = sum(abs(child - blockN*blocks));
numBlocks = blockN;
end
end
numIons = numBlocks*blocks;
else  
end
