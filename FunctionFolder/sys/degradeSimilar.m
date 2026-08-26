function [ranking, bad_rank, fitness] = degradeSimilar(tolerance, ranking1, fitness1, Step)
global POP_STRUC
global ORG_STRUC
bad_rank = 0;
ranking = ranking1;
fitness = fitness1;
it = 1;
while it < length(ranking)
for j = 1 : it-1
if isempty(POP_STRUC.POPULATION(ranking(it)).FINGERPRINT)  || ...
SameStructure(      ranking(j), ranking(it), POP_STRUC) || ...
(POP_STRUC.POPULATION(ranking(it)).Enthalpies(Step) >= 9999)
temp_rank = ranking(it);
for k = it+1:length(ranking)
ranking(k-1) = ranking(k);
end
ranking(end) = temp_rank;
if ORG_STRUC.molecule==1
fitness(ranking(end)) = 100000;
end
bad_rank = bad_rank + 1;
it = it-1;
break;
end
end
if (it + bad_rank == length(ranking))
break;
end
it = it+1;
end
