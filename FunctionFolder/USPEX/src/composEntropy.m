function c_entr = composEntropy()
global POP_STRUC
global ORG_STRUC
N = round((length(POP_STRUC.ranking))*ORG_STRUC.bestFrac);
c_entr = 0;
aver = 0;
entr = 0;
for i = 1:N
for j = i+1:N
f1 = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).numIons;
f1 = f1/sum(f1);
f2 = POP_STRUC.POPULATION(POP_STRUC.ranking(j)).numIons;
f2 = f2/sum(f2);
cos_dist = (dot(f1, f2))/(sqrt(dot(f1, f1)*dot(f2,f2)));
aver = aver + 1 - cos_dist;
end
end
if aver < 0.00001  
c_entr = 0;
else
aver = 1; 
for i = 1:N
for j = i+1:N
f1 = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).numIons;
f1 = f1/sum(f1);
f2 = POP_STRUC.POPULATION(POP_STRUC.ranking(j)).numIons;
f2 = f2/sum(f2);
cos_dist = 1 - (dot(f1, f2))/(sqrt(dot(f1, f1)*dot(f2,f2)));
if abs(1-cos_dist) < 0.001
c_entr = c_entr - (1-cos_dist)/aver;
else
c_entr = c_entr + (1-cos_dist)*log((1-cos_dist)/aver)/aver;
end
end
end
c_entr = -2*c_entr/(N*N-N);
end
