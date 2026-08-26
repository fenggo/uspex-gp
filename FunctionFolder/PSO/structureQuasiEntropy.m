function sQE = structureQuasiEntropy(whichInd, atom_fing)
global POP_STRUC
global ORG_STRUC
sQE = 0;
if ORG_STRUC.molecule==1
numIons=POP_STRUC.POPULATION(whichInd).numMols;
elseif ORG_STRUC.dimension==2
numIons=POP_STRUC.POPULATION(whichInd).Surface_numIons;
else
numIons=POP_STRUC.POPULATION(whichInd).numIons;
end
if (ORG_STRUC.varcomp == 1) || (ORG_STRUC.dimension==2)
weight=ones(1,size(numIons,2));
else
weight = numIons(:)/sum(numIons);
end
[trash, r, c] = size(atom_fing);
for i = 1 : length(numIons)
if numIons(i) > 1
tmp = 0;
k = 0;
start_n = sum(numIons(1:i-1));
for j = start_n+1 : start_n+numIons(i)-1
for j1 = j+1 : start_n+numIons(i)
k = k + 1;
coef1 = 0;
coef2 = 0;
coef3 = 0;
for i1 = 1:r
for i2 = 1:c
coef1 = coef1 + atom_fing(j,i1,i2)*atom_fing(j1,i1,i2)*weight(i1);
coef2 = coef2 + atom_fing(j,i1,i2)*atom_fing(j,i1,i2)*weight(i1);
coef3 = coef3 + atom_fing(j1,i1,i2)*atom_fing(j1,i1,i2)*weight(i1);
end
end
dist = (1-coef1/(sqrt(coef2*coef3)))/2;
if abs(dist-1) < 0.000001
dist = 0.99999;
end
tmp = tmp + (1 - dist)*log(1 - dist);
end
end
sQE = sQE + numIons(i)/sum(numIons)*tmp/k;
end
end
sQE = -sQE;
