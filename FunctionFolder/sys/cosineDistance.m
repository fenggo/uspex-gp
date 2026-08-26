function [dist]=cosineDistance(matrA, matrB, weight)
weight = diag(weight);
coef1 = sum(sum(weight*( matrA.*matrB )));
coef2 = sum(sum(weight*( matrA.*matrA )));
coef3 = sum(sum(weight*( matrB.*matrB )));
dist = (1-coef1/(sqrt(coef2*coef3)))/2;
if abs(dist-1) < 0.000001
dist = 0.99999;
disp(['Fingerprint might be wrong !!!']);
end
