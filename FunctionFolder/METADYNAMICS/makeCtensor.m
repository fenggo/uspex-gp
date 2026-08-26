function [S, C] = makeCtensor()
v = 0.26; 
rules = [1,6,5; 6,2,4; 5,4,3];
C0 = [1-v,v,v,0,0,0; v,1-v,v,0,0,0; v,v,1-v,0,0,0; 0,0,0,0.5-v,0,0; 0,0,0,0,0.5-v,0; 0,0,0,0,0,0.5-v];
S0 = inv(C0);
for i = 1 : 3
for j = 4 : 6
S0(i,j) = S0(i,j)/2;
S0(j,i) = S0(j,i)/2;
end
end
for i = 4 : 6
for j = 4 : 6
S0(i,j) = S0(i,j)/4;
end
end
C = zeros(3,3,3,3);
S = zeros(3,3,3,3);
for i = 1 : 3
for j = 1 : 3
for k = 1 : 3
for l = 1 : 3
C(i,j,k,l) = C0(rules(i,j),rules(k,l));
S(i,j,k,l) = S0(rules(i,j),rules(k,l));
end
end
end
end
tmp = tensorProduct(S,[1 0 0; 0 0 0; 0 0 0]); 
S = S/tmp(1,1);
