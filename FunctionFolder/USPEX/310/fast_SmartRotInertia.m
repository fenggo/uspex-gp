function [Zmatrix, MOLCOORS] = fast_SmartRotInertia(coords, format, num_optFlags, molRot, angleMax, transMax, maxIter)
% FAST_SMARTROTINERTIA — Cython 加速的 SmartRotInertia
%
% 接口与 SmartRotInertia 完全一致:
%   输入: coords (N×3), format (N×3 int), num_optFlags, molRot, angleMax, transMax, maxIter
%   输出: Zmatrix (N×3), MOLCOORS (N×3)
%
% 策略: 将 Z矩阵↔笛卡尔 + 二面角迭代循环交给 Python/Cython,
%       失败时自动回退到原始 SmartRotInertia。

global ORG_STRUC

pythonBin = 'python';
wrapperScript = [ORG_STRUC.USPEXPath '/FunctionFolder/sys/uspex_rotation_wrapper.py'];

% 如果没有 flex_dihedral (即 num_optFlags == 0)，直接用原版
if num_optFlags == 0
    [Zmatrix, MOLCOORS] = SmartRotInertia(coords, format, num_optFlags, molRot, angleMax, transMax, maxIter);
    return;
end

STDMOL = ORG_STRUC.STDMOL;
flex_dihedral = STDMOL(molRot).flex_dihedral;
CN_target = STDMOL(molRot).CN;

% 计算共价半径
N = size(coords, 1);
radiu = zeros(N, 1);
for i = 1:N
    s = covalentRadius(ORG_STRUC.atomType(STDMOL(molRot).types(i)));
    radiu(i) = str2num(s);
end

% 计算惯量主轴
[a, b] = PrincipleAxis(coords);

% --- 准备 JSON 输入 ---
tmpDir = ORG_STRUC.homePath;
tmpId = sprintf('_smrt_%d_%d', round(rand()*1e6), molRot);
tmpInput = [tmpDir '/' tmpId '_input.json'];
tmpMat = [tmpDir '/' tmpId '_result.mat'];

struct_data = struct();
struct_data.coords = coords;
struct_data.format = format;
struct_data.num_optFlags = num_optFlags;
struct_data.flex_dihedral = flex_dihedral;
struct_data.CN_target = CN_target;
struct_data.radii = radiu;
struct_data.inertia_eigvals = diag(b)';
struct_data.inertia_eigvecs = a;
struct_data.angleMax = angleMax;
struct_data.transMax = transMax;
struct_data.maxIter = maxIter;

% 写 JSON
fid = fopen(tmpInput, 'w');
fprintf(fid, '{"structures": [');
fprintf(fid, '{"coords": %s,', mat2json(coords));
fprintf(fid, ' "format": %s,', mat2json(format));
fprintf(fid, ' "num_optFlags": %d,', num_optFlags);
fprintf(fid, ' "flex_dihedral": %s,', mat2json(flex_dihedral));
fprintf(fid, ' "CN_target": %s,', mat2json(CN_target));
fprintf(fid, ' "radii": %s,', mat2json(radiu));
fprintf(fid, ' "inertia_eigvals": %s,', mat2json(diag(b)'));
fprintf(fid, ' "inertia_eigvecs": %s,', mat2json(a));
fprintf(fid, ' "angleMax": %.10f,', angleMax);
fprintf(fid, ' "transMax": %.10f,', transMax);
fprintf(fid, ' "maxIter": %d}', maxIter);
fprintf(fid, ']}');
fclose(fid);

% 调用 Python
cmd = sprintf('%s %s --input=%s --output=%s --nproc=1', ...
              pythonBin, wrapperScript, tmpInput, tmpMat);
[status, msg] = system(cmd);

if status == 0 && exist(tmpMat, 'file')
    try
        tmp = load(tmpMat);
        if tmp.n_structs >= 1
            MOLCOORS = squeeze(tmp.MOLCOORS_all(1, :, :));
            Zmatrix = squeeze(tmp.ZMATRIX_all(1, :, :));
            [nothing, nothing] = system(['rm -f ' tmpDir '/' tmpId '*']);
            return;
        end
    catch
    end
end

% 回退到原版
disp(['fast_SmartRotInertia: Python call failed, falling back to Octave SmartRotInertia']);
[nothing, nothing] = system(['rm -f ' tmpDir '/' tmpId '*']);
[Zmatrix, MOLCOORS] = SmartRotInertia(coords, format, num_optFlags, molRot, angleMax, transMax, maxIter);
end

function s = mat2json(M)
% 将矩阵转为 JSON 数组字符串
s = '[';
if isvector(M) && size(M, 1) == 1
    % 行向量
    for i = 1:length(M)
        if i > 1, s = [s ',']; end
        s = [s sprintf('%.10e', M(i))];
    end
elseif isvector(M)
    % 列向量
    for i = 1:length(M)
        if i > 1, s = [s ',']; end
        s = [s sprintf('%.10e', M(i))];
    end
else
    for i = 1:size(M, 1)
        if i > 1, s = [s ',']; end
        s = [s '['];
        for j = 1:size(M, 2)
            if j > 1, s = [s ',']; end
            s = [s sprintf('%.10e', M(i, j))];
        end
        s = [s ']'];
    end
end
s = [s ']'];
end