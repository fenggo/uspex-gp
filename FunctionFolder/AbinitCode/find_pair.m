function Pair = find_pair(coor, radii)
% find_pair — 矩阵化优化版本
% 原版: O(n²) 双重循环逐对计算距离
% 优化: 一次性计算距离矩阵, 按原始 (i,j) 顺序构建 Pair
%
% 与原版数学等价: 相同的 norm(coor(i,:)-coor(j,:)) < 1.2*(radii(i)+radii(j))
% 相同的 (i,j) 遍历顺序: i=1..n-1, j=i+1..n

N_max = 6;
n_atom = length(radii);

% ---- 矩阵化距离计算 ----
diff = coor - permute(coor, [3 2 1]);
diff = permute(diff, [1 3 2]);
dist = sqrt(sum(diff .^ 2, 3));

% 阈值矩阵
thresh = 1.2 * (radii + radii');

% 全对称连通矩阵
connected_full = dist < thresh;

% ---- 按原始顺序构建 Pair ----
Pair = zeros(n_atom, N_max);

% 遍历上三角 (i,j), i<j, 与原版双重循环顺序一致
for i = 1:n_atom-1
    for j = i+1:n_atom
        if connected_full(i, j)
            cnt_i = Pair(i, N_max);
            if cnt_i < N_max - 1
                Pair(i, cnt_i + 1) = j;
            end
            Pair(i, N_max) = cnt_i + 1;

            cnt_j = Pair(j, N_max);
            if cnt_j < N_max - 1
                Pair(j, cnt_j + 1) = i;
            end
            Pair(j, N_max) = cnt_j + 1;
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
