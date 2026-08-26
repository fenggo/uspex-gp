function findSplit_VC(N, L, R1, R2, splitting)
global ORG_STRUC
S = sum(splitting(1:L)*ORG_STRUC.numIons(1:L,1:end)); 
if N == L+1 
i1 = max(0, ceil((R1-S)/sum(ORG_STRUC.numIons(N,1:end))));
if (S == 0) & (i1 == 0)
i1 = 1;
end
i2 = floor((R2-S)/sum(ORG_STRUC.numIons(N,1:end)));
for i = i1:i2
splitting(N) = i;
ORG_STRUC.splitN = ORG_STRUC.splitN + 1;
if ORG_STRUC.splitN == 1
ORG_STRUC.firstGeneSplit = splitting;
else
ORG_STRUC.firstGeneSplit = vertcat(ORG_STRUC.firstGeneSplit, splitting);
end
end
else
i2 = floor((R2-S)/sum(ORG_STRUC.numIons(L+1,:)));
for i = 0:i2
splitting(L+1) = i;
findSplit_VC(N, L+1, R1, R2, splitting) 
end
end;
ORG_STRUC.firstGeneSplitAll = ORG_STRUC.firstGeneSplit;