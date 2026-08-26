function elasticProperties=calcElasticProperties(elasticMatrix, numIons, atomType, Volume)
h = 6.62606957e-34; 
k = 1.3806488e-23;  
NA=	6.0221413e+23;  
if isempty(elasticMatrix)
elasticProperties=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
return;
end
density = calcDensity(numIons, atomType, Volume);	
C = elasticMatrix;
Kv = ((C(1,1)+C(2,2)+C(3,3)) + 2*(C(1,2)+C(2,3)+C(3,1)))/9;
Gv = ((C(1,1)+C(2,2)+C(3,3)) - (C(1,2)+C(2,3)+C(3,1)) + 3*(C(4,4)+C(5,5)+C(6,6)))/15;
Ev = 1/((1/(3*Gv))+(1/(9*Kv)));
vv  = 0.5*(1-((3*Gv)/(3*Kv+Gv)));
S = inv(C);      
Kr = 1/((S(1,1)+S(2,2)+S(3,3))+2*(S(1,2)+S(2,3)+S(3,1)));
Gr = 15/(4*(S(1,1)+S(2,2)+S(3,3))-4*(S(1,2)+S(2,3)+S(3,1))+3*(S(4,4)+S(5,5)+S(6,6)));
Er = 1/((1/(3*Gr))+(1/(9*Kr)));
vr = 0.5*(1-((3*Gr)/(3*Kr+Gr)));
K_h = (Kv+Kr)/2;   
G_h = (Gv+Gr)/2;   
E_h = (Ev+Er)/2;   
v_h = (vv+vr)/2;   
GK  = G_h/K_h;     
E = eig(C) ; 
is_stable = all(E>0);  
if G_h*K_h < 1e-3
Hv = NaN;
else
Hv = 2*((G_h*G_h*G_h)/(K_h*K_h))^0.585-3;   
end
V_p = ((K_h+(4*G_h/3))/density)^0.5;         
V_s = (G_h/density)^0.5;                     
V_m = (1/3*(2/(V_s^3)+1/(V_p^3)))^(-1/3);    
atomicMass=zeros(1, length(atomType));
for i = 1:length(atomType)
atomicMass(i) = elementMass(atomType(i));
end	  
n = sum(numIons);                               
Sum_M = sum(numIons.*atomicMass);                   
De_T = h/k*((3*n)/(4*pi)*(NA*density/Sum_M))^(1/3)*V_m*100000;    
Kg =0 ;
elasticProperties=[K_h, G_h, E_h, v_h, GK, Hv, Kg, De_T, V_m, V_s, V_p, is_stable];
