function fitness = FitnessRankingCorrection(fitness)
global POP_STRUC
global ORG_STRUC
[nothing, ranking] = sort(fitness);
POP_STRUC.ranking = ranking;
bad_rank = 0;
Step = length(ORG_STRUC.abinitioCode);
[ranking, bad_rank, fitness] = degradeSimilar(ORG_STRUC.toleranceFing, ranking, fitness, Step);
POP_STRUC.fitness = fitness;
POP_STRUC.ranking = ranking;
POP_STRUC.bad_rank = bad_rank;
save ('Current_POP.mat', 'POP_STRUC')
