function [order, FINGERPRINT, atom_fing, V, success, soap_fp] = fast_fingerprint_fullatom(LATTICE, COORDINATES, numIons, atomType, Intra_map)
% FAST_FINGERPRINT_FULLATOM - 全原子指纹（带分子内距离过滤，Cython 加速）
%
% 对应 ReadJobs_310.m 行58-64 的全原子指纹计算:
%   makeMatrices + Intra_MOL_dist 过滤 + fingerprint_calc
%
% 用法:
%   [order, FINGERPRINT, atom_fing, V, success] = fast_fingerprint_fullatom(LATTICE, COORDINATES, numIons, atomType, Intra_map)
%   [order, FINGERPRINT, atom_fing, V, success, soap_fp] = fast_fingerprint_fullatom(..., Intra_map)
%
% Intra_map: N×N 矩阵，0=分子内（过滤掉），1=分子间（保留）
% soap_fp:   SOAP 指纹向量（1D），若计算失败则为空

global ORG_STRUC

success = 0;
order = []; FINGERPRINT = []; atom_fing = []; V = 0;
soap_fp = [];

Rmax = 12.0; sigma = 0.05; delta = 0.08;
if isfield(ORG_STRUC, 'RmaxFing'), Rmax = ORG_STRUC.RmaxFing; end
if isfield(ORG_STRUC, 'sigmaFing'), sigma = ORG_STRUC.sigmaFing; end
if isfield(ORG_STRUC, 'deltaFing'), delta = ORG_STRUC.deltaFing; end

% SOAP 参数（可被 ORG_STRUC.soapFing = 1 启用）
useSoap = false;
if isfield(ORG_STRUC, 'soapFing') && ORG_STRUC.soapFing
    useSoap = true;
end
soapRCut = 6.0; soapNMax = 8; soapLMax = 6;
if isfield(ORG_STRUC, 'soapRCut'), soapRCut = ORG_STRUC.soapRCut; end
if isfield(ORG_STRUC, 'soapNMax'), soapNMax = ORG_STRUC.soapNMax; end
if isfield(ORG_STRUC, 'soapLMax'), soapLMax = ORG_STRUC.soapLMax; end

pythonBin = 'python';
wrapperScript = [ORG_STRUC.USPEXPath '/FunctionFolder/sys/uspex_fingerprint_wrapper.py'];

tmpDir = ORG_STRUC.homePath;
tmpId = sprintf('_fp_%d', round(rand()*1e6));
tmpPoscar = [tmpDir '/' tmpId '_poscar'];
tmpMat = [tmpDir '/' tmpId '_result.mat'];
tmpIntraMat = [tmpDir '/' tmpId '_intra.mat'];

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

% 写 Intra_map 为 .mat 文件（-mat7-binary 保证 scipy 可读）
save('-mat7-binary', tmpIntraMat, 'Intra_map');

% 构建 Python 命令
soapFlag = '';
if useSoap
  soapFlag = sprintf(' --soap --soap-r-cut=%f --soap-n-max=%d --soap-l-max=%d', soapRCut, soapNMax, soapLMax);
end
cmd = sprintf('%s %s --poscar=%s --output=%s --rmax=%f --sigma=%f --delta=%f --intra-map=%s%s', ...
              pythonBin, wrapperScript, tmpPoscar, tmpMat, Rmax, sigma, delta, tmpIntraMat, soapFlag);
[status, msg] = system(cmd);

if status ~= 0
  disp('fast_fingerprint_fullatom: Python failed, falling back');
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
  if isfield(tmp, 'soap_fp')
    soap_fp = tmp.soap_fp;
  end
else
  disp('fast_fingerprint_fullatom: .mat file not found');
end

% 清理
[nothing, nothing] = system(['rm -f ' tmpDir '/' tmpId '*']);
