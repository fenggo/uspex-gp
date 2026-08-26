function [new_Coord, N_Moved]= move_atom_Mutation(Ind_No)
global ORG_STRUC
global POP_STRUC
new_Lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
new_Coord = POP_STRUC.POPULATION(Ind_No).COORDINATES;
if length(new_Lattice)==6
new_Lattice = latConverter(new_Lattice);
end
temp_potLat = latConverter(new_Lattice);
tmp = zeros(sum(POP_STRUC.POPULATION(Ind_No).numIons),1);
for i = 1 : sum(POP_STRUC.POPULATION(Ind_No).numIons)
tmp(i) = i;
end
atoms = sort_order(POP_STRUC.POPULATION(Ind_No).order, tmp, 'ascend');
howManyMut = round(ORG_STRUC.howManyMut*sum(POP_STRUC.POPULATION(Ind_No).numIons));
if howManyMut < 1
howManyMut = 1;
end
N_Moved = RandInt(1,1,[1,howManyMut]);
for i = 1 : N_Moved
deviat_dist = randn(3,1)*2;
new_Coord(atoms(i),1) = new_Coord(atoms(i),1) + deviat_dist(1)/temp_potLat(1);
new_Coord(atoms(i),2) = new_Coord(atoms(i),2) + deviat_dist(2)/temp_potLat(2);
new_Coord(atoms(i),3) = new_Coord(atoms(i),3) + deviat_dist(3)/temp_potLat(3);
new_Coord(atoms(i),1) = new_Coord(atoms(i),1) - floor(new_Coord(atoms(i),1));
new_Coord(atoms(i),2) = new_Coord(atoms(i),2) - floor(new_Coord(atoms(i),2));
new_Coord(atoms(i),3) = new_Coord(atoms(i),3) - floor(new_Coord(atoms(i),3));
end
