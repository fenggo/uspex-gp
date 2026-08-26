function [new_Coord] = move_all_mol_Mutation(MOLECULES, numMols, LATTICE, order)
% move_all_mol_Mutation — 向量化版本
% 原版: 逐分子循环 randn(3,1) → 位移 → wrap
% 优化: 一次性 randn(N,3) → 广播位移 → 一次性 wrap
% 与原版数学等价: 相同的 koef 公式、相同的 randn 分布、相同的 floor wrap

global ORG_STRUC
max_sigma = ORG_STRUC.howManyMut;
N = sum(ORG_STRUC.numMols);

% 收集所有分子质心 (Cartesian)
new_Coord = zeros(N, 3);
for ind = 1:N
    new_Coord(ind, :) = MOLECULES(ind).ZMATRIX(1, :);
end

% 转为分数坐标
new_Lattice = LATTICE;
if length(new_Lattice) == 6
    new_Lattice = latConverter(new_Lattice);
end
new_Coord = new_Coord / new_Lattice;
temp_potLat = latConverter(new_Lattice);

% 按 order 排序
[~, ranking] = sort(order);
r1 = order(ranking(1));
rN = order(ranking(N));

if rN > r1
    % ---- 向量化: 一次性计算所有 koef ----
    rI_all = order(ranking);                     % N×1
    koef_all = (rN - rI_all) / (rN - r1);        % N×1

    % ---- 一次性生成所有随机位移 (N×3) ----
    deviat_dist = randn(N, 3) .* (max_sigma * koef_all);

    % ---- 转为分数坐标位移 (广播除法) ----
    deviat_frac = deviat_dist ./ temp_potLat(1:3);

    % ---- 一次性应用位移 ----
    new_Coord(ranking, :) = new_Coord(ranking, :) + deviat_frac;

    % ---- 一次性 wrap 到 [0,1) ----
    new_Coord(ranking, :) = new_Coord(ranking, :) - floor(new_Coord(ranking, :));
end
end
