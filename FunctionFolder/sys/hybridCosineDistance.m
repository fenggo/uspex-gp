function [dist] = hybridCosineDistance(fpA, fpB, soapA, soapB, weight, alpha)
% HYBRIDCOSINEDISTANCE - 混合 RDF + SOAP 余弦距离
%
% 用法:
%   dist = hybridCosineDistance(fpA, fpB, soapA, soapB, weight, alpha)
%
% 输入:
%   fpA, fpB   - RDF 指纹矩阵 (n_species^2 x n_bins)，同 cosineDistance
%   soapA, soapB - SOAP 指纹向量 (1D)，可为空 []
%   weight     - RDF 权重向量 (同 cosineDistance 的 weight)
%   alpha      - RDF 权重 [0,1]，默认 0.9
%
% 返回:
%   dist       - 混合余弦距离 [0, 1]
%                0 = 完全相同, 1 = 完全不同
%
% 当 soapA 或 soapB 为空时，退化为纯 RDF cosineDistance。
%
% 原理: RDF 捕获长程堆积，SOAP 捕获局部角度环境。
%       两者 Spearman 相关仅 ~0.24，高度互补。
%       混合距离 alpha*RDF + (1-alpha)*SOAP 能区分
%       RDF 相似但分子取向不同的多形性结构。

if nargin < 6
    alpha = 0.9;
end

% 如果 SOAP 指纹缺失，退化为纯 RDF 距离
if isempty(soapA) || isempty(soapB)
    dist = cosineDistance(fpA, fpB, weight);
    return;
end

% --- RDF 余弦距离 ---
w = diag(weight);
coef1 = sum(sum(w * (fpA .* fpB)));
coef2 = sum(sum(w * (fpA .* fpA)));
coef3 = sum(sum(w * (fpB .* fpB)));
if coef2 < 1e-30 || coef3 < 1e-30
    rdf_dist = 1.0;
else
    rdf_cos = coef1 / sqrt(coef2 * coef3);
    rdf_cos = max(min(rdf_cos, 1.0), -1.0);
    rdf_dist = (1 - rdf_cos) / 2;  % 归一化到 [0, 1]
end

% --- SOAP 余弦距离 ---
sa = soapA(:);
sb = soapB(:);
na = sqrt(sa' * sa);
nb = sqrt(sb' * sb);
if na < 1e-30 || nb < 1e-30
    soap_dist = 1.0;
else
    soap_cos = (sa' * sb) / (na * nb);
    soap_cos = max(min(soap_cos, 1.0), -1.0);
    soap_dist = (1 - soap_cos) / 2;  % 归一化到 [0, 1]
end

% --- 混合距离 ---
dist = alpha * rdf_dist + (1 - alpha) * soap_dist;
