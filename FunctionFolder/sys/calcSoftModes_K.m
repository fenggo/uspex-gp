function [freq, eigvector] = calcSoftModes_K(R_val, N_val, val, lat, coords, numIons, kVector0)
global ORG_STRUC
rec_lat = zeros(3,3);
rec_lat(1,:) = 2*pi*cross(lat(2,:), lat(3,:))/det(lat);
rec_lat(2,:) = 2*pi*cross(lat(3,:), lat(1,:))/det(lat);
rec_lat(3,:) = 2*pi*cross(lat(1,:), lat(2,:))/det(lat);
kVector = kVector0*rec_lat;
tic
checkConnectivity = 1;
N_i = sum(numIons);
D = zeros(3*N_i, 3*N_i);
vect = zeros(1,3);
bonds_add = zeros(1,7);
same_bond = 0.05; 
max_bond = 5; 
small_bond = -0.37*log(ORG_STRUC.goodBonds); 
nL = 2; 
if length(lat) == 6
lat = latConverter(lat);
end
at_types = zeros(1,N_i);
for k = 1 : N_i
tmp = k;
while tmp > 0
at_types(k) = at_types(k) + 1;
tmp = tmp - numIons(at_types(k));
end
end
N_bonds1 = 0; 
colors = zeros(2*nL+1, 2*nL+1, 2*nL+1, N_i); 
Ncolors = ((2*nL)^3)* N_i;
for k1 = 1 : 2*nL+1
for k2 = 1 : 2*nL+1
for k3 = 1 : 2*nL+1
for i = 1 : N_i
colors(k1,k2,k3,i) = i + ((2*nL+1)^2)*N_i*(k1-1) + (2*nL+1)*N_i*(k2-1) + N_i*(k3-1);
end
end
end
end
colors1 = colors;
it = 1;
for i = 1 : N_i    
for j = i : N_i
for k1 = -1 : 1
for k2 = -1 : 1
for k3 = -1 : 1
if (i == j) & (k1 == 0) & (k2 == 0) & (k3 == 0)
continue;
end
vect(1) = coords(i,1) + k1 - coords(j,1); 
vect(2) = coords(i,2) + k2 - coords(j,2);
vect(3) = coords(i,3) + k3 - coords(j,3);
delta = sqrt(sum((vect*lat).^2)) - R_val(at_types(i)) - R_val(at_types(j)); 
if delta < max_bond
N_bonds1 = N_bonds1 + 1;
bonds_add(1) = i;
bonds_add(2) = j;
bonds_add(3) = delta;
bonds_add(4) = 0; 
bonds_add(5) = k1;
bonds_add(6) = k2;
bonds_add(7) = k3;
if N_bonds1 == 1
bonds = bonds_add;
else
bonds = vertcat(bonds, bonds_add);
for bb = 1 : N_bonds1-1
if bonds(bb,3) > bonds_add(3)
for bbb = N_bonds1 : -1 : bb+1
bonds(bbb,:) = bonds(bbb-1,:);  
end
bonds(bb,:) = bonds_add;
break;
end
end
end
end
end 
end  
end   
end    
end     
bond_types = 0;
N_bonds = 0;
for bb = 1 : N_bonds1
if bonds(bb,4) == 0;
bond_types = bond_types + 1;
bonds(bb,4) = bond_types;
N_bonds = N_bonds + 1;
if N_bonds == 1
bonds_sorted = bonds(bb,:);
else
bonds_sorted = vertcat(bonds_sorted, bonds(bb,:));       
end
for bbb = bb+1 : N_bonds1
if bonds(bbb,3) > bonds(bb,3) + same_bond
break;
end
if (bonds(bbb,4) == 0) & (at_types(bonds(bbb,1)) == at_types(bonds(bb,1))) & (at_types(bonds(bbb,2)) == at_types(bonds(bb,2)))
bonds(bbb,4) = bond_types;
N_bonds = N_bonds + 1;
bonds_sorted = vertcat(bonds_sorted, bonds(bbb,:));      
end
end
end
end
bondTypesDensity = zeros(1,bond_types);
h_tmp = ones(1,bond_types);
it = 0;
for bt = 1 : bond_types
while (bonds_sorted(it+1,4) == bt) & (it < N_bonds)
it = it + 1;
a = at_types(bonds_sorted(it,1));
b = at_types(bonds_sorted(it,2));
bond_used = 0;
for m1 = 1 : 2*nL+1
for m2 = 1 : 2*nL+1
for m3 = 1 : 2*nL+1
n1 = m1 + bonds_sorted(it,5);
n2 = m2 + bonds_sorted(it,6);
n3 = m3 + bonds_sorted(it,7);
if (n1 > 2*nL+1) | (n1 < 1) | (n2 > 2*nL+1) | (n2 < 1) | (n3 > 2*nL+1) | (n3 < 1)
continue;
end
if (colors(nL+1,nL+1,nL+1,bonds_sorted(it,2)) == colors(nL+1+bonds_sorted(it,5), nL+1+bonds_sorted(it,6), nL+1+bonds_sorted(it,7),bonds_sorted(it,1))) & (bonds_sorted(it,3) > small_bond(a,b))
continue;
end
if (checkConnectivity == 0) & (bonds_sorted(it,3) > small_bond(a,b))
continue;
end
bond_used = 1;
cj = colors1(m1,m2,m3,bonds_sorted(it,2));
ci = colors1(n1,n2,n3,bonds_sorted(it,1));
if checkConnectivity == 1
for r1 = 1 : 2*nL+1
for r2 = 1 : 2*nL+1
for r3 = 1 : 2*nL+1
for m = 1 : N_i
if colors1(r1,r2,r3,m) == cj
colors1(r1,r2,r3,m) = ci;
end
end
end
end
end
end
end
end
end
bondTypesDensity(bt) = bondTypesDensity(bt) + bond_used;
if bond_used == 0
bonds_sorted(it,1) = 0;
bonds_sorted(it,2) = 0;    
end
if it == N_bonds
break;
end
end
if checkConnectivity == 1
for i = 1 : N_i
for j = i+1 : N_i
if colors1(nL+1,nL+1,nL+1,i) == colors1(nL+1,nL+1,nL+1,j)   
for m1 = 1 : 2*nL+1
for m2 = 1 : 2*nL+1
for m3 = 1 : 2*nL+1
if colors1(m1,m2,m3,i) ~= colors1(m1,m2,m3,j)
cj = colors1(m1,m2,m3,i);
ci = colors1(m1,m2,m3,j);
for r1 = 1 : 2*nL+1
for r2 = 1 : 2*nL+1
for r3 = 1 : 2*nL+1
for m = 1 : N_i
if colors1(r1,r2,r3,m) == cj
colors1(r1,r2,r3,m) = ci;
end
end
end
end
end
end
end
end
end
end
end
end
end 
if it+1 < N_bonds
if bonds_sorted(it+1,3) < max(max(small_bond))
continue;
end
end
colors = colors1;
connected = 1;
c1 = colors(1,1,1,1);
c2 = 0;
for k1 = nL : nL+2
for k2 = nL : nL+2
for k3 = nL : nL+2
for i = 1 : N_i
if (nL == 0) & (colors(1,1,1,i) ~= c1)
connected = 0;
c2 = 1;
else      
if (colors(k1,k2,k3,i) ~= c1) & (c2 == 0)
c2 = colors(k1,k2,k3,i);
elseif (colors(k1,k2,k3,i) ~= c1) & (colors(k1,k2,k3,i) ~= c2)
connected = 0;
end
end
end
end
end
end
if connected | (it >= N_bonds) | (checkConnectivity == 0)
for dumm = it+1 : N_bonds
bonds_sorted(dumm,1) = 0;
bonds_sorted(dumm,2) = 0;
end
break;
end
end
nu_factor = zeros(1,N_i);
for k = 1 : N_i
nu_full = 0;
for m = 1 : N_bonds
if (bonds_sorted(m,1) == k) | (bonds_sorted(m,2) == k)
nu_full = nu_full + exp(-1*bonds_sorted(m,3)/0.37);
end
end
nu_factor(k) = val(at_types(k))/nu_full;
end
it = 0;
for bt = 1 : bond_types
while (bonds_sorted(it+1,4) == bt) & (it < N_bonds)
it = it + 1;
if bonds_sorted(it,1) == 0
if it == N_bonds
break;
end
continue;
end
a = at_types(bonds_sorted(it,1));
b = at_types(bonds_sorted(it,2));
delta = bonds_sorted(it,3);
R_a = R_val(a) + delta/2;
R_b = R_val(b) + delta/2;
if (R_a < 0.05)
R_b = R_b - 0.05 + R_a;
R_a = 0.05; 
end
if (R_b < 0.05)
R_a = R_a - 0.05 + R_b;
R_b = 0.05;
end
nu = exp(-delta/0.37);
EN_a = 0.481*N_val(a)/R_a;
EN_b = 0.481*N_val(b)/R_b;
CN_a = val(a)/(nu*nu_factor(bonds_sorted(it,1)));
CN_b = val(b)/(nu*nu_factor(bonds_sorted(it,2)));
f_ab = 0.25*abs(EN_a - EN_b)/sqrt(EN_a*EN_b);
X_ab = sqrt((EN_a*EN_b)/(CN_a*CN_b));
b = bonds_sorted(it,1); 
a = bonds_sorted(it,2);
k1 = bonds_sorted(it,5);
k2 = bonds_sorted(it,6);
k3 = bonds_sorted(it,7);
vect(1) = coords(b,1) + k1 - coords(a,1); 
vect(2) = coords(b,2) + k2 - coords(a,2);
vect(3) = coords(b,3) + k3 - coords(a,3);
dR_ab1 = vect*lat;
dR_ab2 = -1*dR_ab1;
cos_x = (dR_ab1(1))/norm(dR_ab1);
cos_y = (dR_ab1(2))/norm(dR_ab1);
cos_z = (dR_ab1(3))/norm(dR_ab1);
cos_x = round(cos_x*1000000)/1000000;
cos_y = round(cos_y*1000000)/1000000;
cos_z = round(cos_z*1000000)/1000000;
i = sqrt(-1);
phase_k1 = exp(i*dot(kVector, dR_ab1));
phase_k2 = exp(i*dot(kVector, dR_ab2));
H = X_ab*exp(-2.7*f_ab);
if a == b
if ~((k1 == 0) & (k2 == 0) & (k3 == 0)) 
D((a-1)*3+1, (a-1)*3+1) = D((a-1)*3+1, (a-1)*3+1) + H*(1 - phase_k1)*cos_x*cos_x; 
D((a-1)*3+1, (a-1)*3+2) = D((a-1)*3+1, (a-1)*3+2) + H*(1 - phase_k1)*cos_x*cos_y; 
D((a-1)*3+1, (a-1)*3+3) = D((a-1)*3+1, (a-1)*3+3) + H*(1 - phase_k1)*cos_x*cos_z; 
D((a-1)*3+2, (a-1)*3+1) = D((a-1)*3+2, (a-1)*3+1) + H*(1 - phase_k1)*cos_y*cos_x; 
D((a-1)*3+2, (a-1)*3+2) = D((a-1)*3+2, (a-1)*3+2) + H*(1 - phase_k1)*cos_y*cos_y; 
D((a-1)*3+2, (a-1)*3+3) = D((a-1)*3+2, (a-1)*3+3) + H*(1 - phase_k1)*cos_y*cos_z; 
D((a-1)*3+3, (a-1)*3+1) = D((a-1)*3+3, (a-1)*3+1) + H*(1 - phase_k1)*cos_z*cos_x; 
D((a-1)*3+3, (a-1)*3+2) = D((a-1)*3+3, (a-1)*3+2) + H*(1 - phase_k1)*cos_z*cos_y; 
D((a-1)*3+3, (a-1)*3+3) = D((a-1)*3+3, (a-1)*3+3) + H*(1 - phase_k1)*cos_z*cos_z; 
end
else 
D((a-1)*3+1, (a-1)*3+1) = D((a-1)*3+1, (a-1)*3+1) + H*cos_x*cos_x; 
D((a-1)*3+1, (a-1)*3+2) = D((a-1)*3+1, (a-1)*3+2) + H*cos_x*cos_y; 
D((a-1)*3+1, (a-1)*3+3) = D((a-1)*3+1, (a-1)*3+3) + H*cos_x*cos_z; 
D((a-1)*3+1, (b-1)*3+1) = D((a-1)*3+1, (b-1)*3+1) - H*phase_k1*cos_x*cos_x; 
D((a-1)*3+1, (b-1)*3+2) = D((a-1)*3+1, (b-1)*3+2) - H*phase_k1*cos_x*cos_y; 
D((a-1)*3+1, (b-1)*3+3) = D((a-1)*3+1, (b-1)*3+3) - H*phase_k1*cos_x*cos_z; 
D((a-1)*3+2, (a-1)*3+1) = D((a-1)*3+2, (a-1)*3+1) + H*cos_y*cos_x; 
D((a-1)*3+2, (a-1)*3+2) = D((a-1)*3+2, (a-1)*3+2) + H*cos_y*cos_y; 
D((a-1)*3+2, (a-1)*3+3) = D((a-1)*3+2, (a-1)*3+3) + H*cos_y*cos_z; 
D((a-1)*3+2, (b-1)*3+1) = D((a-1)*3+2, (b-1)*3+1) - H*phase_k1*cos_y*cos_x; 
D((a-1)*3+2, (b-1)*3+2) = D((a-1)*3+2, (b-1)*3+2) - H*phase_k1*cos_y*cos_y; 
D((a-1)*3+2, (b-1)*3+3) = D((a-1)*3+2, (b-1)*3+3) - H*phase_k1*cos_y*cos_z; 
D((a-1)*3+3, (a-1)*3+1) = D((a-1)*3+3, (a-1)*3+1) + H*cos_z*cos_x; 
D((a-1)*3+3, (a-1)*3+2) = D((a-1)*3+3, (a-1)*3+2) + H*cos_z*cos_y; 
D((a-1)*3+3, (a-1)*3+3) = D((a-1)*3+3, (a-1)*3+3) + H*cos_z*cos_z; 
D((a-1)*3+3, (b-1)*3+1) = D((a-1)*3+3, (b-1)*3+1) - H*phase_k1*cos_z*cos_x; 
D((a-1)*3+3, (b-1)*3+2) = D((a-1)*3+3, (b-1)*3+2) - H*phase_k1*cos_z*cos_y; 
D((a-1)*3+3, (b-1)*3+3) = D((a-1)*3+3, (b-1)*3+3) - H*phase_k1*cos_z*cos_z; 
D((b-1)*3+1, (b-1)*3+1) = D((b-1)*3+1, (b-1)*3+1) + H*cos_x*cos_x; 
D((b-1)*3+1, (b-1)*3+2) = D((b-1)*3+1, (b-1)*3+2) + H*cos_x*cos_y; 
D((b-1)*3+1, (b-1)*3+3) = D((b-1)*3+1, (b-1)*3+3) + H*cos_x*cos_z; 
D((b-1)*3+1, (a-1)*3+1) = D((b-1)*3+1, (a-1)*3+1) - H*phase_k2*cos_x*cos_x; 
D((b-1)*3+1, (a-1)*3+2) = D((b-1)*3+1, (a-1)*3+2) - H*phase_k2*cos_x*cos_y; 
D((b-1)*3+1, (a-1)*3+3) = D((b-1)*3+1, (a-1)*3+3) - H*phase_k2*cos_x*cos_z; 
D((b-1)*3+2, (b-1)*3+1) = D((b-1)*3+2, (b-1)*3+1) + H*cos_y*cos_x; 
D((b-1)*3+2, (b-1)*3+2) = D((b-1)*3+2, (b-1)*3+2) + H*cos_y*cos_y; 
D((b-1)*3+2, (b-1)*3+3) = D((b-1)*3+2, (b-1)*3+3) + H*cos_y*cos_z; 
D((b-1)*3+2, (a-1)*3+1) = D((b-1)*3+2, (a-1)*3+1) - H*phase_k2*cos_y*cos_x; 
D((b-1)*3+2, (a-1)*3+2) = D((b-1)*3+2, (a-1)*3+2) - H*phase_k2*cos_y*cos_y; 
D((b-1)*3+2, (a-1)*3+3) = D((b-1)*3+2, (a-1)*3+3) - H*phase_k2*cos_y*cos_z; 
D((b-1)*3+3, (b-1)*3+1) = D((b-1)*3+3, (b-1)*3+1) + H*cos_z*cos_x; 
D((b-1)*3+3, (b-1)*3+2) = D((b-1)*3+3, (b-1)*3+2) + H*cos_z*cos_y; 
D((b-1)*3+3, (b-1)*3+3) = D((b-1)*3+3, (b-1)*3+3) + H*cos_z*cos_z; 
D((b-1)*3+3, (a-1)*3+1) = D((b-1)*3+3, (a-1)*3+1) - H*phase_k2*cos_z*cos_x; 
D((b-1)*3+3, (a-1)*3+2) = D((b-1)*3+3, (a-1)*3+2) - H*phase_k2*cos_z*cos_y; 
D((b-1)*3+3, (a-1)*3+3) = D((b-1)*3+3, (a-1)*3+3) - H*phase_k2*cos_z*cos_z; 
end 
if it == N_bonds
break;
end
end
end
[eigvector, freq] = eig(D);
freq = real(freq);
