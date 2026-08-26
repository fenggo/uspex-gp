function [order, FINGERPRINT, atom_fing, V, success] = fast_fingerprint(LATTICE, COORDINATES, numIons, atomType)
% FAST_FINGERPRINT - 分子质心指纹（Cython 加速）
%
% 对应 ReadJobs_310.m 行51-52 的质心指纹计算:
%   makeMatrices + fingerprint_calc
%
% 用法:
%   [order, FINGERPRINT, atom_fing, V, success] = fast_fingerprint(LATTICE, COORDINATES, numIons, atomType)
%
% 不传 Intra_map, 不做分子内距离过滤

global ORG_STRUC

success = 0;
order = []; FINGERPRINT = []; atom_fing = []; V = 0;

Rmax = 12.0; sigma = 0.05; delta = 0.08;
if isfield(ORG_STRUC, 'RmaxFing'), Rmax = ORG_STRUC.RmaxFing; end
if isfield(ORG_STRUC, 'sigmaFing'), sigma = ORG_STRUC.sigmaFing; end
if isfield(ORG_STRUC, 'deltaFing'), delta = ORG_STRUC.deltaFing; end

pythonBin = 'python';
wrapperScript = [ORG_STRUC.USPEXPath '/FunctionFolder/sys/uspex_fingerprint_wrapper.py'];

tmpDir = ORG_STRUC.homePath;
tmpId = sprintf('_fp_%d', round(rand()*1e6));
tmpPoscar = [tmpDir '/' tmpId '_poscar'];
tmpMat = [tmpDir '/' tmpId '_result.mat'];

% 写 POSCAR
fid = fopen(tmpPoscar, 'w');
fprintf(fid, 'tmp\n1.0000\n');
for i = 1:3
  fprintf(fid, '  %14.10f  %14.10f  %14.10f\n', LATTICE(i,1), LATTICE(i,2), LATTICE(i,3));
end
atomSyms = '';
for i = 1:length(atomType)
  atomSyms = [atomSyms megaDoof(atomType(i)) ' '];
end
fprintf(fid, '%s\n', atomSyms);
numIonsStr = '';
for i = 1:length(numIons)
  numIonsStr = [numIonsStr num2str(numIons(i)) ' '];
end
fprintf(fid, '%s\nDirect\n', numIonsStr);
for i = 1:size(COORDINATES, 1)
  fprintf(fid, '  %14.10f  %14.10f  %14.10f\n', COORDINATES(i,1), COORDINATES(i,2), COORDINATES(i,3));
end
fclose(fid);

% 调用 Python (不传 --intra-map, 不做分子内过滤)
cmd = sprintf('%s %s --poscar=%s --output=%s --rmax=%f --sigma=%f --delta=%f', ...
              pythonBin, wrapperScript, tmpPoscar, tmpMat, Rmax, sigma, delta);
[status, msg] = system(cmd);

if status ~= 0
  disp('fast_fingerprint: Python failed, falling back');
  disp(msg);
  [nothing, nothing] = system(['rm -f ' tmpDir '/' tmpId '*']);
  return;
end

% load .mat
if exist(tmpMat, 'file')
  tmp = load(tmpMat);
  order = tmp.order;
  FINGERPRINT = tmp.FINGERPRINT;
  atom_fing = tmp.atom_fing;
  V = tmp.V;
  success = 1;
else
  disp('fast_fingerprint: .mat file not found');
end

% 清理
[nothing, nothing] = system(['rm -f ' tmpDir '/' tmpId '*']);
