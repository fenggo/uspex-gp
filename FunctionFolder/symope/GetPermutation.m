function [permutation, permutationBack] = GetPermutation(nsym);
if (nsym > 2) & (nsym < 16) 
if rand > 0.5
permutation = [3 2 1]; 
permutationBack = permutation;
else
permutation = [1 2 3]; 
permutationBack = permutation;
end
elseif (nsym > 15) & (nsym < 75) 
ListPerm = [1 2 3; 3 2 1; 1 3 2; 2 1 3; 3 1 2; 2 3 1];  
id = RandInt(1,1,6)+1;
permutation = ListPerm(id,:); 
if id == 5
permutationBack = ListPerm(6,:);  
elseif id == 6
permutationBack = ListPerm(5,:);  
else
permutationBack = permutation;
end
else
permutation = [1 2 3];   
permutationBack = permutation;
end
