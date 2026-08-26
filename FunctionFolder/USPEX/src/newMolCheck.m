function goodDist = newMolCheck(MOLECULES, LATTICE, MtypeLIST, miniMatrix)
% newMolCheck — 矩阵化优化版本
% 将逐对循环替换为矩阵运算，以空间换时间，不改变计算逻辑。
%
% 优化策略:
%   Part 1 (分子间): 预计算所有原子对的分数坐标差矩阵 →
%                   最小镜像 → 笛卡尔距离矩阵 → 阈值矩阵一次性比较
%   Part 2 (分子内): 预计算 124 偏移矩阵 → 广播计算所有原子对到
%                   所有周期镜像的距离 → 逐原子对取 min → 阈值比较
%
% 与原版数学等价: 相同的 inv(LATTICE)/round/norm 路径。

global ORG_STRUC;
STDMOL = ORG_STRUC.STDMOL;
nMols = length(MtypeLIST);

% ===================================================================
% 预计算: 收集所有分子的原子坐标和类型
% ===================================================================
molCoords = cell(1, nMols);
molTypes  = cell(1, nMols);
molNAtoms = zeros(1, nMols);

for i = 1:nMols
    molCoords{i} = MOLECULES(i).MOLCOORS;
    molTypes{i}  = STDMOL(MtypeLIST(i)).types(:);
    molNAtoms(i) = length(molTypes{i});
end

invLat = inv(LATTICE);

% ===================================================================
% Part 1: 分子间距离检查 (矩阵化)
% 原版: 四重循环 + 逐对 CalcDist
% 优化: 每对分子一次性计算距离矩阵 + 阈值矩阵
% ===================================================================
for i = 1:nMols - 1
    coords_i = molCoords{i};
    types_i  = molTypes{i};
    n_i = molNAtoms(i);
    if n_i == 0, continue; end

    frac_i = coords_i * invLat;  % n_i × 3

    for k = i+1:nMols
        coords_k = molCoords{k};
        types_k  = molTypes{k};
        n_k = molNAtoms(k);
        if n_k == 0, continue; end

        frac_k = coords_k * invLat;  % n_k × 3

        % 所有原子对的分数坐标差: n_i × n_k × 3
        diff_frac = frac_i - permute(frac_k, [3 2 1]);
        diff_frac = permute(diff_frac, [1 3 2]);

        % 最小镜像约定
        diff_frac = diff_frac - round(diff_frac);

        % 转为笛卡尔坐标差并计算欧氏距离
        % diff_cart(a,b,:) = squeeze(diff_frac(a,b,:))' * LATTICE
        diff_frac_2d = reshape(diff_frac, n_i * n_k, 3);
        cart_diff = diff_frac_2d * LATTICE;
        cart_diff = reshape(cart_diff, n_i, n_k, 3);
        dists = sqrt(sum(cart_diff .^ 2, 3));  % n_i × n_k

        % 阈值矩阵: thresh(a,b) = miniMatrix(type_i(a), type_k(b))
        thresh = miniMatrix(types_i, types_k);  % n_i × n_k

        % 一次性比较
        if any(any(dists <= thresh))
            goodDist = 0;
            return;
        end
    end
end

% ===================================================================
% Part 2: 分子内距离检查 (矩阵化)
% 原版: 对每个原子对调用 bsxfun + pdist2 计算 124 个镜像距离
% 优化: 预计算 124 偏移矩阵, 广播计算所有原子对 × 所有镜像距离
% ===================================================================

% 预计算 124 个周期偏移向量 (5×5×5 - 原点)
[X1, Y1, Z1] = deal(-2:2);
[X2, Y2, Z2] = meshgrid(X1, Y1, Z1);
offsetMat = [X2(:), Y2(:), Z2(:)];          % 125 × 3
offsetMat(all(offsetMat == 0, 2), :) = [];  % 124 × 3
offsetCart = offsetMat * LATTICE;            % 124 × 3

for ind = 1:nMols
    coords = molCoords{ind};
    types  = molTypes{ind};
    nAtoms = molNAtoms(ind);

    if nAtoms < 2
        continue;
    end

    for n = 2:nAtoms
        % 原子 n 的所有 124 个周期镜像: 124 × 3
        images_n = coords(n, :) + offsetCart;

        % 原子 1..n-1: (n-1) × 3
        atoms_m = coords(1:n-1, :);

        % 广播计算所有 (m, 镜像) 对的距离
        % diff: (n-1) × 124 × 3
        diff = reshape(atoms_m, n-1, 1, 3) - reshape(images_n, 1, 124, 3);
        dists = sqrt(sum(diff .^ 2, 3));  % (n-1) × 124

        % 每个 m 取最小距离 (最近镜像)
        minDists = min(dists, [], 2);  % (n-1) × 1

        % 阈值向量: thresh(m) = miniMatrix(type(m), type(n))
        thresh = miniMatrix(types(1:n-1), types(n));

        % 一次性比较
        if any(minDists <= thresh)
            goodDist = 0;
            return;
        end
    end
end

goodDist = 1;
end
