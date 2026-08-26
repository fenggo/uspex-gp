function [freq, eigvector] = fast_calcSoftModes(N_val, val, Ind_No, POP_STRUC)
% FAST_CALCSOFTMODES — Cython 加速的 calcSoftModes_molecules
%
% 接口与 calcSoftModes_molecules 完全一致:
%   输入: N_val (价电子数), val (价数), Ind_No, POP_STRUC
%   输出: freq (3N×3N 对角矩阵), eigvector (3N×3N)
%
% 策略: 将键搜索+动力学矩阵+eig 的完整计算交给 Python/Cython,
%       Octave 只负责数据准备和结果回传。
%       失败时自动回退到原始 calcSoftModes_molecules。

global ORG_STRUC

pythonBin = 'python';
wrapperScript = [ORG_STRUC.USPEXPath '/FunctionFolder/sys/uspex_softmode_wrapper.py'];

% --- 从 POP_STRUC 提取结构数据 (与 calcSoftModes_molecules 完全一致) ---
N_i = sum(POP_STRUC.POPULATION(Ind_No).numIons);
lat = POP_STRUC.POPULATION(Ind_No).LATTICE;
coords = POP_STRUC.POPULATION(Ind_No).COORDINATES;

if length(lat) == 6
    lat = latConverter(lat);
end

% 原子类型 (1-based -> 0-based for Python)
at_types = zeros(N_i, 1);
for k = 1:N_i
    tmp = k;
    while tmp > 0
        at_types(k) = at_types(k) + 1;
        tmp = tmp - POP_STRUC.POPULATION(Ind_No).numIons(at_types(k));
    end
end
at_types = at_types - 1;  % 0-based

% 共价半径
R_val = zeros(length(ORG_STRUC.atomType), 1);
for i = 1:length(ORG_STRUC.atomType)
    s = covalentRadius(ceil(ORG_STRUC.atomType(i)));
    R_val(i) = str2num(s);
end

% goodBonds
goodBonds = ORG_STRUC.goodBonds;

% val_arr 传入 (确保为列向量)
val_arr = val(:);
N_val_arr = N_val(:);

% --- 写临时 .mat 文件 ---
tmpDir = ORG_STRUC.homePath;
tmpId = sprintf('_csm_%d', round(rand()*1e6));
tmpInput  = [tmpDir '/' tmpId '_input.mat'];
tmpOutput = [tmpDir '/' tmpId '_output.mat'];

save('-mat7-binary', tmpInput, 'coords', 'lat', 'at_types', ...
     'R_val', 'N_val_arr', 'val_arr', 'goodBonds');

% --- 调用 Python ---
cmd = sprintf('%s %s --input=%s --output=%s', pythonBin, wrapperScript, tmpInput, tmpOutput);
[status, msg] = system(cmd);

% 清理输入文件
[nothing, nothing] = system(['rm -f ' tmpInput]);

if status ~= 0
    % Python 失败: 回退到原始函数
    disp('fast_calcSoftModes: Python failed, falling back to calcSoftModes_molecules');
    disp(msg);
    [nothing, nothing] = system(['rm -f ' tmpOutput]);
    [freq, eigvector] = calcSoftModes_molecules(N_val, val, Ind_No, POP_STRUC);
    return;
end

% --- 读取结果 ---
if exist(tmpOutput, 'file')
    result = load(tmpOutput);
    [nothing, nothing] = system(['rm -f ' tmpOutput]);

    freq = result.freq;           % (3N × 3N) 对角矩阵
    eigvector = result.eigvector; % (3N × 3N)

    % 验证维度
    if size(freq, 1) ~= 3*N_i || size(eigvector, 1) ~= 3*N_i
        disp('fast_calcSoftModes: dimension mismatch, falling back');
        [freq, eigvector] = calcSoftModes_molecules(N_val, val, Ind_No, POP_STRUC);
    end

    % fprintf('fast_calcSoftModes: %.3fs\n', result.time_used);
else
    disp('fast_calcSoftModes: output file not found, falling back');
    [freq, eigvector] = calcSoftModes_molecules(N_val, val, Ind_No, POP_STRUC);
end
