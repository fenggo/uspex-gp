function [decomposable, samecomp, fitness] = CheckDecomposition(convex_hull, E_hull, A, Energy)
warning off
N_Comp        = size(convex_hull, 1);
N_Block       = size(convex_hull, 2);
decomposable  = 1;
samecomp      = 0;
fitness = [];
  
disp('check compositions ');
for i = 1:N_Comp
    if sameComposition(A, convex_hull(i,:))
       fitness = Energy - E_hull(i);
       samecomp = 1;
       break;
    end
end

if isempty(fitness)
    Comp_List     = nchoosek([1:N_Comp],N_Block);
    N_Comp_List   = nchoosek(N_Comp,N_Block);
    C = zeros(N_Block, N_Block);
    E = zeros(N_Block, 1);
    form_Eng = [];
    N_paths = 0;
    for i = 1: N_Comp_List
        for j = 1: N_Block
            C(j,:) = convex_hull(Comp_List(i,j), :);
            C(j,:) = C(j,:)/sum(C(j,:));  
            E(j)   = E_hull(Comp_List(i,j));
        end
        if abs(det(C)) > 1.0e-8  
            X = A/C;
            if isempty(find(X<0)) 
            N_paths = N_paths + 1;
            form_Eng(N_paths) = Energy - X*E/sum(A);  
            end
        end
    end

    if N_paths == 0
        disp(['Odd compositions: ' num2str(A), ' unable to find the decomposition paths']);
        fitness = 10000;
    else
        fitness = max(form_Eng);
    end
end

if fitness < 0
   decomposable = 0;
end
