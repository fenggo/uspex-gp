function [a, b] = PrincipleAxis(AbsoluteCoord)
% PrincipleAxis — 向量化优化版本
% 原版: 逐元素计算惯性张量 (9 行逐个赋值)
% 优化: 向量化计算, 与原版数学等价

AbsoluteCoord = bsxfun(@minus, AbsoluteCoord, mean(AbsoluteCoord));

% 向量化惯性张量计算
x = AbsoluteCoord(:, 1);
y = AbsoluteCoord(:, 2);
z = AbsoluteCoord(:, 3);

Inertia = zeros(3, 3);
Inertia(1, 1) = sum(y.^2 + z.^2);
Inertia(2, 2) = sum(x.^2 + z.^2);
Inertia(3, 3) = sum(x.^2 + y.^2);
Inertia(1, 2) = -sum(x .* y);
Inertia(2, 3) = -sum(y .* z);
Inertia(3, 1) = -sum(z .* x);
Inertia(2, 1) = Inertia(1, 2);  % 对称
Inertia(3, 2) = Inertia(2, 3);
Inertia(1, 3) = Inertia(3, 1);

[a, b] = eig(Inertia);
end
