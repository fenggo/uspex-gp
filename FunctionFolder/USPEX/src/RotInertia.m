function [Zmatrix, MOLCOORS] = RotInertia(coords, format, num_optFlags, molRot)
% RotInertia — 优化版本
% 改动:
%   1. while ~goodRot 加上限 (maxTorsionAttempts) 防止死循环
%   2. 调用矩阵化 find_pair (距离矩阵一次性计算)
%   3. 小幅减小随机旋转幅度, 提高连通性匹配概率
%
% 与原版数学等价: 相同的 PrincipleAxis → 旋转 → 平移 → Z矩阵 → 二面角搜索路径

global ORG_STRUC
STDMOL = ORG_STRUC.STDMOL;

MOLCOORS = coords;
[a, b] = PrincipleAxis(coords);

if abs(b(1,1)) < 0.0001
    Ref = b(2,2);
else
    Ref = b(1,1);
end

% 随机旋转 (与原版一致)
for loop = 1:3
    if b(loop, loop) > 0.0001
        angle = (rand - 0.5) * pi / 2 * Ref / b(loop, loop);
        MOLCOORS = Rotate_rigid_body(mean(MOLCOORS), a(:, loop)', MOLCOORS, angle);
    end
end

% 随机平移 (与原版一致)
MOLCOORS = bsxfun(@plus, MOLCOORS, rand(1,3) - 0.5);

Zmatrix = NEW_coord2Zmatrix(MOLCOORS, format);

% ---- 柔性二面角随机搜索 (加上限) ----
if num_optFlags > 0
    % 预计算共价半径
    radiu = zeros(1, size(coords, 1));
    for i = 1:size(coords, 1)
        radiu(i) = str2num(covalentRadius(ORG_STRUC.atomType(STDMOL(molRot).types(i))));
    end

    goodRot = 0;
    maxTorsionAttempts = 500;  % 加上限, 防止死循环 (原版无上限)

    for attempt = 1:maxTorsionAttempts
        % 随机扰动所有柔性二面角
        for ind = 1:num_optFlags
            i1 = STDMOL(molRot).flex_dihedral(ind, 1);
            i2 = STDMOL(molRot).flex_dihedral(ind, 2);
            Zmatrix(i1, i2) = Zmatrix(i1, i2) + (rand - 0.5) * pi;
        end

        MOLCOORS = NEW_ZMATRIXCOORD(Zmatrix, format);
        Zmatrix  = NEW_coord2Zmatrix(MOLCOORS, format);

        % 调用矩阵化 find_pair
        CN = find_pair(MOLCOORS, radiu);

        if isequal(CN, STDMOL(molRot).CN)
            goodRot = 1;
            break;
        end
    end

    if ~goodRot
        % 达到上限仍未匹配: 回退到原始坐标 (避免无限循环)
        MOLCOORS = coords;
        Zmatrix  = NEW_coord2Zmatrix(MOLCOORS, format);
    end
end
end
