function [fitness, convex_hull0] = phase_update_hull(N_Block, N_Point, numBlocks, Energy, numIons)
tic
fitness = zeros(1, N_Point);
convex_hull = 1000000*ones(N_Block, N_Block + 2);
ch_add      = zeros(1, N_Block + 2);
finished    = zeros(1, N_Point);
for i = 1 : N_Point
for j = 1: N_Block
if numBlocks(i,j) == sum(numBlocks(i,:))
finished(i) = 1;
if Energy(i) <= convex_hull(j, N_Block+1) + 0.0000001
convex_hull(j, 1 : N_Block)   = numBlocks(i,:);  
convex_hull(j, N_Block + 1)   = Energy(i);       
convex_hull(j, N_Block + 2)   = i;               
end
end
end
end
for i = 1 : N_Point
if finished(i)
continue;
end
[decomposable, samecomp, fitness(i)] = CheckDecomposition(convex_hull(:, 1:N_Block), ...
convex_hull(:, N_Block+1), numBlocks(i,:), Energy(i));
if decomposable == 0 
ch_add  = [numBlocks(i,:) Energy(i) i];
if samecomp == 1  
for j = 1:size(convex_hull, 1)
if sameComposition(numBlocks(i,:), convex_hull(j,:))
convex_hull(j,:) = ch_add;
break
end
end
else              
convex_hull   = [convex_hull; ch_add];
end
convex_hull1 = convex_hull;
convex_hull(N_Block+1:end, :) = [];
for j = N_Block+1 : size(convex_hull1, 1)
tmp_hull = convex_hull1;
tmp_hull(j,:) = [];
[tmp1, tmp2] = CheckDecomposition(tmp_hull(:, 1:N_Block), ...
tmp_hull(:, N_Block+1), convex_hull1(j, 1:N_Block), convex_hull1(j, 1+N_Block));
if tmp1 == 0  
convex_hull = [convex_hull; convex_hull1(j,:)];
end
end
end
end
for i = 1 : N_Point
[decomposable, samecomp, fitness(i)] = CheckDecomposition(convex_hull(:, 1:N_Block), ...
convex_hull(:, N_Block+1), numBlocks(i,:), Energy(i));
end
convex_hull0 = Transform_Convex_Hull(convex_hull, numIons, 0);
t = toc;
%disp(['Calculating fitness in ' num2str(t) ' sec'])
%disp(' ')
