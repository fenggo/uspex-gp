function update_PSO_GENERATION(IND, fitness, type)
global PSO_STRUC
global POP_STRUC
PSO_STRUC.GENERATION(POP_STRUC.generation).quasiEntropy = oldQuasiEntropy();
PSO_STRUC.GENERATION(POP_STRUC.generation).Fitness = min(fitness);
for i=1:length(IND)
ID = POP_STRUC.POPULATION(IND(i)).Number;
PSO_STRUC.GENERATION(POP_STRUC.generation).BestID(i) = ID;
end
[nothing, ranking] = sort(POP_STRUC.DoneOrder);
for i = 1 : length(fitness)
if POP_STRUC.DoneOrder(ranking(i)) > 0
ID = POP_STRUC.POPULATION(ranking(i)).Number;
PSO_STRUC.POPULATION(ID).Fitness = fitness(ranking(i));
end
end
