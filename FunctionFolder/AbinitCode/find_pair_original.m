function Pair = find_pair(coor, radii)
% find_pair — 矩阵化优化版本
% 原版: O(n²) 双重循环逐对计算距离
% 优化: 一次性计算距离矩阵 + 阈值矩阵, 空间换时间
%
% 与原版数学等价: 相同的 norm(coor(i,:)-coor(j,:)) < 1.2*(radii(i)+radii(j))

N_max = 6;
n_atom = length(radii);

% ---- 矩阵化距离计算 ----
% diff(i,j,:) = coor(i,:) - coor(j,:),  n_atom × n_atom × 3
diff = coor - permute(coor, [3 2 1]);
diff = permute(diff, [1 3 2]);

% dist(i,j) = norm(diff(i,j,:))
dist = sqrt(sum(diff .^ 2, 3));

% 阈值矩阵: thresh(i,j) = 1.2 * (radii(i) + radii(j))
thresh = 1.2 * (radii + radii');

% 连通性矩阵: connected(i,j) = dist(i,j) < thresh(i,j) 且 i < j
connected = (dist < thresh) & triu(true(n_atom), 1);

% ---- 构建 Pair 矩阵 (与原版格式完全一致) ----
Pair = zeros(n_atom, N_max);

for i = 1:n_atom
    % 找到与原子 i 相连的所有原子 j
    js = find(connected(i, :));
    for idx = 1:length(js)
        j = js(idx);
        Pair(i, N_max) = Pair(i, N_max) + 1;
        if Pair(i, N_max) <= N_max - 1
            Pair(i, Pair(i, N_max)) = j;
        end
    end
end

% 孤立原子检查 (与原版一致)
if n_atom > 1
    for i = 1:n_atom
        if (Pair(i, N_max) == 0) && (N_max > 1)
            disp(['atom_' num2str(i) ' is not connected to any other atom']);
            disp(['Please check your MOL file again. Serious WARNING.... ']);
        end
    end
end
end
