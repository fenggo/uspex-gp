% test_newMolCheck.m — 验证矩阵化 newMolCheck 与原版数学等价性
%
% 用法: 在 USPEX 310 计算目录下运行
%   matlab -nodisplay -r "test_newMolCheck; exit"

clear; clc;

% 加载 USPEX 全局变量
global ORG_STRUC;
load('ORG_STRUC.mat');

fprintf('=== newMolCheck 矩阵化版本验证 ===\n');
fprintf('atomType: %s\n', strjoin(ORG_STRUC.atomType, ' '));

% 读取 INPUT.txt 中的 miniMatrix
fid = fopen('INPUT.txt', 'r');
miniMatrix = [];
while ~feof(fid)
    line = fgetl(fid);
    if contains(line, '% minDistMatIce')
        miniMatrix = zeros(4,4);
        for i = 1:4
            line = fgetl(fid);
            miniMatrix(i,:) = sscanf(line, '%f %f %f %f');
        end
        break;
    end
end
fclose(fid);

if isempty(miniMatrix)
    % 使用默认值
    miniMatrix = [1.74 1.69 1.71 1.29;
                  1.69 1.64 1.66 1.24;
                  1.71 1.66 1.69 1.26;
                  1.29 1.24 1.26 1.02];
end
fprintf('miniMatrix:\n'); disp(miniMatrix);

% 生成测试数据: 随机分子坐标
rng(42);
nMols = 8;
LATTICE = [25 0 0; 0 25 0; 0 0 25];
MtypeLIST = [1 1 1 1 2 2 2 2];

% 创建模拟 MOLECULES 结构
MOLECULES = struct('MOLCOORS', cell(1, nMols));
STDMOL(1).types = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1];  % TNT-like: 21 atoms
STDMOL(2).types = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];  % CL20-like: 36 atoms
ORG_STRUC.STDMOL = STDMOL;

for i = 1:nMols
    nAtoms = length(STDMOL(MtypeLIST(i)).types);
    MOLECULES(i).MOLCOORS = rand(nAtoms, 3) * 20;  % 随机坐标
end

% 测试 100 组随机配置
fprintf('\n测试 100 组随机分子配置...\n');
nTests = 100;
agree = 0;
disagree = 0;
total_time_old = 0;
total_time_new = 0;

for t = 1:nTests
    % 重新生成随机坐标
    for i = 1:nMols
        nAtoms = length(STDMOL(MtypeLIST(i)).types);
        MOLECULES(i).MOLCOORS = rand(nAtoms, 3) * 20;
    end

    % 运行原版 (需要原始 newMolCheck_original.m)
    tic;
    r1 = newMolCheck_original(MOLECULES, LATTICE, MtypeLIST, miniMatrix);
    t_old = toc;
    total_time_old = total_time_old + t_old;

    % 运行新版
    tic;
    r2 = newMolCheck(MOLECULES, LATTICE, MtypeLIST, miniMatrix);
    t_new = toc;
    total_time_new = total_time_new + t_new;

    if r1 == r2
        agree = agree + 1;
    else
        disagree = disagree + 1;
        fprintf('  ✗ 测试 %d: 原版=%d, 新版=%d\n', t, r1, r2);
    end
end

fprintf('\n=== 结果 ===\n');
fprintf('一致: %d/%d (%.1f%%)\n', agree, nTests, agree/nTests*100);
fprintf('不一致: %d\n', disagree);
fprintf('原版总耗时: %.3f s (平均 %.3f ms)\n', total_time_old, total_time_old/nTests*1000);
fprintf('新版总耗时: %.3f s (平均 %.3f ms)\n', total_time_new, total_time_new/nTests*1000);
if total_time_old > 0
    fprintf('加速比: %.1fx\n', total_time_old / total_time_new);
end

% 额外: 大规模压力测试
fprintf('\n=== 大规模压力测试 (1000 次) ===\n');
nLarge = 1000;
t_start = tic;
for t = 1:nLarge
    for i = 1:nMols
        nAtoms = length(STDMOL(MtypeLIST(i)).types);
        MOLECULES(i).MOLCOORS = rand(nAtoms, 3) * 20;
    end
    r = newMolCheck(MOLECULES, LATTICE, MtypeLIST, miniMatrix);
end
t_large = toc(t_start);
fprintf('1000 次新版耗时: %.3f s (平均 %.3f ms)\n', t_large, t_large/nLarge*1000);
