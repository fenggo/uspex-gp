function update_USPEX_GENERATION(IND, fitness, type)
global USPEX_STRUC
global POP_STRUC
USPEX_STRUC.GENERATION(POP_STRUC.generation).quasiEntropy = oldQuasiEntropy();
USPEX_STRUC.GENERATION(POP_STRUC.generation).Fitness = min(fitness);
if type == 2  
USPEX_STRUC.GENERATION(POP_STRUC.generation).convex_hull= POP_STRUC.convex_hull;
USPEX_STRUC.GENERATION(POP_STRUC.generation).composEntropy= composEntropy();
end
ID = [];
for i=1:length(IND)
if IND(i) > 0
ID = [ID  POP_STRUC.POPULATION(IND(i)).Number];
end
end
USPEX_STRUC.GENERATION(POP_STRUC.generation).BestID = ID;
USPEX_STRUC.GENERATION(POP_STRUC.generation).ID = length(USPEX_STRUC.POPULATION); 
[nothing, ranking] = sort(POP_STRUC.DoneOrder);
for i = 1 : length(fitness)
if POP_STRUC.DoneOrder(ranking(i)) > 0
ID = POP_STRUC.POPULATION(ranking(i)).Number;
USPEX_STRUC.POPULATION(ID).Fitness = fitness(ranking(i));
end
end
