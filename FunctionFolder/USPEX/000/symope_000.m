function [candidate, newLattice, errorS] = symope_000(nsym, numIonsFull, lattice, minDistMatrix)
ellipse_mode = 1;
E = [1.0 0 0; 0 1 0; 0 0 1];                                    
I = [-1.0, 0.0, 0.0; 0.0, -1.0, 0.0; 0.0, 0.0, -1.0];       
C2x = [-1 0 0; 0 1 0; 0 0 1]; 
C2y = [1 0 0; 0 -1 0; 0 0 1]; 
C2z = [1 0 0; 0 1 0; 0 0 -1]; 
Hz = [1 0 0; 0 1 0; 0 0 -1];  
errorS = 0;
newLattice = lattice;
volLat = det(lattice);
numIons = sum(numIonsFull);
nI = numIons;
minDistance = minDistMatrix(1,1); 
if strcmp(nsym, 'E')
candidate = rand(sum(numIons),3) - 0.5;
elseif strcmp(nsym, 'Ci') | strcmp(nsym, 'S2')  
candidate = [];
while 1
if numIons == 1;   
tmp = [0 0 0];          
candidate = cat(1,candidate,tmp); 
break;
end
if numIons == 0 
break
end
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5; 
end
if sum(abs(tmp)) < 0.01
continue
end;
nextAtom = tmp*I;
cand = cat(1,tmp,nextAtom);
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
numIons = numIons - 2;
end
elseif ~isempty(findstr(nsym, 'I')) | ~isempty(findstr(nsym, 'i'))   
newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
candidate = [];
if mod(numIons,2) == 1 
candidate = [0 0 0];
numIons = numIons - 1;
end
if (numIons < 60) & (mod(numIons,12) ~= 0) & (numIons ~= 20) & (numIons ~= 30) & (numIons ~= 32)
if (numIons < 40) | (numIons == 58) | (numIons == 46)
status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
errorS = 1;
%    unix(['echo ' status ' > error_cluster_symmetry']);
end
end
n60 = 0;
numIons1 = numIons;
while 1
numIons1 = numIons1 - 60;
if (numIons1 ~= 0) & (numIons1 < 60)
if (numIons1 ~= 12) & (numIons1 ~= 20) & (numIons1 ~= 24) & (numIons1 ~= 30) & (numIons1 ~= 32) & (numIons1 ~= 36) 
if (numIons1 < 40) | (numIons1 == 58) | (numIons1 == 46)
break;
end
end
else
n60 = n60 + 1;
end
end
GR = (sqrt(5)+1)/2;  
alpha = atan(1/GR); 
R1 = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
R5 = [1 0 0 ; 0 cos(2*pi/5) -sin(2*pi/5); 0 sin(2*pi/5) cos(2*pi/5)]; 
n = 1;
while n <= n60  
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5; 
end
if (strcmp(nsym, 'Ih') | strcmp(nsym, 'ih')) 
if (mod(n60,2) == 1) & (n == n60)  
tmp(3) = 0;   
else 
tmp = cat(1,tmp,tmp*[1 0 0; 0 1 0; 0 0 -1]);
end
end       
cand = tmp;
for i = 1 : 4
nextAtom = tmp*inv(R1)*(R5^i)*R1; 
cand = cat(1,cand,nextAtom);
end
nextAtoms1 = cand*[0 1 0; 0 0 1; 1 0 0]; 
nextAtoms2 = cand*[0 0 1; 1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtoms1);
cand = cat(1,cand,nextAtoms2);
nextAtoms = cand*[-1 0 0; 0 -1 0; 0 0 1]; 
cand = cat(1,cand,nextAtoms);
nextAtoms = cand*[1 0 0; 0 -1 0; 0 0 -1]; 
cand = cat(1,cand,nextAtoms);
if (strcmp(nsym, 'Ih') | strcmp(nsym, 'ih')) & ~((mod(n60,2) == 1) & (n == n60))  
numIons = numIons - 120;
n = n + 2;
else
numIons = numIons - 60;
n = n + 1;
end        
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if n > n60 
break
end
end
if numIons > 0  
for x = 0 : 10
for y = 0 : 6
for z = 0 : 4
if 12*x + 20*y + 30*z == numIons
break
end
end
if 12*x + 20*y + 30*z == numIons
break
end
end
if 12*x + 20*y + 30*z == numIons
break
end
end
for i = 1 : x 
tmp = [1 -1/GR 0]*rand*0.5;
while ellipse_mode & (norm(tmp) > 0.5)
tmp = [1 -1/GR 0]*rand*0.5;
end
cand = tmp;
cand = cat(1,cand,tmp*[0 1 0; 0 0 1; 1 0 0]); 
cand = cat(1,cand,tmp*[0 0 1; 1 0 0; 0 1 0]); 
nextAtoms = cand*[-1 0 0; 0 -1 0; 0 0 1]; 
cand = cat(1,cand,nextAtoms);
nextAtoms = cand*[1 0 0; 0 -1 0; 0 0 -1]; 
cand = cat(1,cand,nextAtoms);        
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end     
end
for i = 1 : y 
tmp1 = rand*0.5;
while abs(tmp1) < 0.05
tmp1 = rand - 0.5;
end
tmp = [tmp1 -tmp1 tmp1];
while ellipse_mode & (norm(tmp) > 0.5)
tmp1 = rand*0.5;
while abs(tmp1) < 0.05
tmp1 = rand - 0.5;
end
tmp = [tmp1 -tmp1 tmp1];
end
cand = tmp;
nextAtoms = cand*[-1 0 0; 0 -1 0; 0 0 1]; 
cand = cat(1,cand,nextAtoms);
nextAtoms = cand*[-1 0 0; 0 1 0; 0 0 -1]; 
cand = cat(1,cand,nextAtoms);
candTMP = cand;
for i = 1 : 4
nextAtom = candTMP*inv(R1)*(R5^i)*R1; 
cand = cat(1,cand,nextAtom);
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end        
end
for i = 1 : z 
tmp1 = rand - 0.5;
while abs(tmp1) < 0.05
tmp1 = rand - 0.5;
end
tmp = [0 0 tmp1]; 
cand = tmp;
nextAtoms = cand*[1 0 0; 0 -1 0; 0 0 -1]; 
cand = cat(1,cand,nextAtoms); 
nextAtoms1 = cand*[0 1 0; 0 0 1; 1 0 0]; 
nextAtoms2 = cand*[0 0 1; 1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtoms1);
cand = cat(1,cand,nextAtoms2);
candTMP = cand;
for i = 1 : 4
nextAtom = candTMP*inv(R1)*(R5^i)*R1; 
cand = cat(1,cand,nextAtom);
end     
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end 
end
end
elseif ~isempty(findstr(nsym, 'S')) | ~isempty(findstr(nsym, 's')) 
newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
candidate = [];
rotRank = str2num(nsym(2:end));   
while 1
if numIons == 1;   
tmp = [0 0 0];          
candidate = cat(1,candidate,tmp); 
break;
end
if numIons < rotRank        
cand = [0 0 rand-0.5];
for i = 2 : numIons
tmp = [0 0 rand-0.5]; 
cand = cat(1,cand,tmp); 
end 
numIons = 0;
else
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5; 
end
cand = tmp;
tooClose = 0;      
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat*(Hz^i);
if norm((nextAtom-tmp)*newLattice) < minDistance
tooClose = 1;
end
cand = cat(1,cand,nextAtom);
end
if tooClose
numIons = numIons - 1;
cand = [0 0 tmp(3)];
else
numIons = numIons - rotRank;
end
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if numIons == 0 
break
end
end
elseif strcmp(nsym, 'Th') | strcmp(nsym, 'th')   
newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
candidate = [];
if (numIons ~= 1) & (numIons ~= 6) & (numIons ~= 7) & (numIons ~= 8) & (numIons ~= 9) & (numIons < 12)
status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
[nothing, nothing] = unix(['echo ' status ' > error_cluster_symmetry']);
errorS = 1;
end
if mod(numIons,2) == 1 
candidate = [0 0 0];
numIons = numIons - 1;
end
n12 = 0;
numIons1 = numIons;
while 1
numIons1 = numIons1 - 12;
if (numIons1 ~= 6) & (numIons1 ~= 8) & (numIons1 ~= 0) & (numIons1 < 12)
break;
else
n12 = n12 + 1;
end
end
n = 1;
while n <= n12  
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if (mod(n12,2) == 1) & (n == n12) 
tmp(ceil(rand*3)) = 0;  
end   
cand = tmp;
nextAtom = tmp*[-1 0 0; 0 -1 0; 0 0 1]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[1 0 0; 0 -1 0; 0 0 -1];  
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[-1 0 0; 0 1 0; 0 0 -1];  
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 1 0; 0 0 1; 1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 1; 1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 1 0; 0 0 -1; -1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 -1; 1 0 0; 0 -1 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 -1 0; 0 0 1; -1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 -1; -1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtom);     
nextAtom = tmp*[0 -1 0; 0 0 -1; 1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 1; -1 0 0; 0 -1 0]; 
cand = cat(1,cand,nextAtom); 
if ~((mod(n12,2) == 1) & (n == n12))
numIons = numIons - 24;
n = n + 2;
cand1 = -1*cand; 
cand = cat(1,cand,cand1);
else
numIons = numIons - 12;
n = n + 1;
end     
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if n > n12 
break
end
end
if numIons > 0  
n1 = ceil(numIons/8);
n2 = floor(numIons/6);
n = n1 + round(rand*(n2-n1)); 
x = 4*n - round(numIons/2);
y = round(numIons/2) - 3*n;
for i = 1 : x   
tmp = rand - 0.5;
while abs(tmp) < 0.05
tmp = rand - 0.5;
end
candidate = cat(1, candidate, [0 0 tmp]);
candidate = cat(1, candidate, [0 0 -1*tmp]);
candidate = cat(1, candidate, [0 tmp 0]);
candidate = cat(1, candidate, [0 -1*tmp 0]);
candidate = cat(1, candidate, [tmp 0 0]);
candidate = cat(1, candidate, [-1*tmp 0 0]);
end
for i = 1 : y  
tmp = rand - 0.5;
while abs(tmp) < 0.05 | (ellipse_mode & (norm([tmp tmp tmp]) > 0.5))
tmp = rand - 0.5;
end
candidate = cat(1, candidate, [tmp tmp tmp]);
candidate = cat(1, candidate, [tmp tmp -1*tmp]);
candidate = cat(1, candidate, [tmp -1*tmp tmp]);
candidate = cat(1, candidate, [tmp -1*tmp -1*tmp]);
candidate = cat(1, candidate, [-1*tmp tmp tmp]);
candidate = cat(1, candidate, [-1*tmp tmp -1*tmp]);
candidate = cat(1, candidate, [-1*tmp -1*tmp tmp]);
candidate = cat(1, candidate, [-1*tmp -1*tmp -1*tmp]);
end
end    
elseif ~isempty(findstr(nsym, 'T')) | ~isempty(findstr(nsym, 't'))   
newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
candidate = [];
if (numIons ~= 1) & (numIons < 4)
status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
[nothing, nothing] = unix(['echo ' status ' > error_cluster_symmetry']);
errorS = 1;
end
if mod(numIons,2) == 1 
candidate = [0 0 0];
numIons = numIons - 1;
end
n12 = 0;
numIons1 = numIons;
while 1
numIons1 = numIons1 - 12;
if (numIons1 ~= 0) & (numIons1 < 4)
break;
else
n12 = n12 + 1;
end
end
n = 1;
while n <= n12  
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if (strcmp(nsym, 'Td') | strcmp(nsym, 'td')) & (mod(n12,2) == 1) & (n == n12)  
tmp(1) = tmp(2);   
end   
cand = tmp;
nextAtom = tmp*[-1 0 0; 0 -1 0; 0 0 1]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[1 0 0; 0 -1 0; 0 0 -1];  
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[-1 0 0; 0 1 0; 0 0 -1];  
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 1 0; 0 0 1; 1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 1; 1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 1 0; 0 0 -1; -1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 -1; 1 0 0; 0 -1 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 -1 0; 0 0 1; -1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 -1; -1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtom);     
nextAtom = tmp*[0 -1 0; 0 0 -1; 1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 1; -1 0 0; 0 -1 0]; 
cand = cat(1,cand,nextAtom); 
if (strcmp(nsym, 'Td') | strcmp(nsym, 'td')) & ~((mod(n12,2) == 1) & (n == n12))  
numIons = numIons - 24;
n = n + 2;
cand1 = cand; 
cand1(:,1) = -1*cand(:,2);
cand1(:,2) = -1*cand(:,1);
cand = cat(1,cand,cand1);
else
numIons = numIons - 12;
n = n + 1;
end     
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if n > n12 
break
end
end
if numIons > 0  
n1 = ceil(numIons/6);
n2 = floor(numIons/4);
n = n1 + round(rand*(n2-n1)); 
x = 3*n - round(numIons/2);
y = round(numIons/2) - 2*n;
for i = 1 : x   
tmp = rand - 0.5;
while abs(tmp) < 0.05  | (ellipse_mode & (norm([tmp tmp tmp]) > 0.5))
tmp = rand - 0.5;
end
candidate = cat(1, candidate, [tmp tmp tmp]);
candidate = cat(1, candidate, [tmp -1*tmp -1*tmp]);
candidate = cat(1, candidate, [-1*tmp -1*tmp tmp]);
candidate = cat(1, candidate, [-1*tmp tmp -1*tmp]);
end
for i = 1 : y  
tmp = rand - 0.5;
while abs(tmp) < 0.05
tmp = rand - 0.5;
end
candidate = cat(1, candidate, [0 0 tmp]);
candidate = cat(1, candidate, [0 0 -1*tmp]);
candidate = cat(1, candidate, [0 tmp 0]);
candidate = cat(1, candidate, [0 -1*tmp 0]);
candidate = cat(1, candidate, [tmp 0 0]);
candidate = cat(1, candidate, [-1*tmp 0 0]);
end
end    
elseif ~isempty(findstr(nsym, 'Cv')) | ~isempty(findstr(nsym, 'cv')) 
rotRank = str2num(nsym(3:end));   
if rotRank > 2
newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
end
candidate = [];
while 1
if numIons == 1;   
tmp = [0 0 0];          
candidate = cat(1,candidate,tmp); 
break;
end
if numIons < 2*rotRank 
if numIons < rotRank      
cand = [0 0 rand-0.5];
for i = 2 : numIons
tmp = [0 0 rand-0.5]; 
cand = cat(1,cand,tmp); 
end
numIons = 0;
else                   
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
r = norm(tmp);
tmp = [sqrt(r^2-tmp(3)^2) 0 tmp(3)];  
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
numIons = numIons - rotRank;
end
else
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
cand = tmp;
if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  
numIons = numIons - 1;
cand = [0 0 tmp(3)];          
else
tooClose = 0;
for i = 1 : rotRank        
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
if norm((nextAtom-tmp)*newLattice) < minDistance 
tooClose = 1;
tmp = (tmp + nextAtom)/2;
cand = tmp;
break;
end
end
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
for i = 1 : rotRank
if tooClose
break;
end
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end 
if tooClose
numIons = numIons - rotRank;
else
numIons = numIons - 2*rotRank;
end
end
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if numIons == 0 
break
end
end
elseif ~isempty(findstr(nsym, 'Ch')) | ~isempty(findstr(nsym, 'ch')) 
rotRank = str2num(nsym(3:end));   
if rotRank > 2
newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
end
candidate = [];
while 1
if numIons == 1;   
tmp = [0 0 0];          
candidate = cat(1,candidate,tmp); 
break;
end
if numIons < 2*rotRank 
if numIons < rotRank      
cand = [0 0 rand-0.5];
while abs(cand(3)) < 0.01
cand = [0 0 rand-0.5]; 
end 
cand = cat(1,cand,cand*Hz);    
numIons = numIons - 2;
else                   
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm([tmp(1) tmp(2)]) > 0.5)
tmp = rand(1,3) - 0.5;
end
tmp(3) = 0;
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
numIons = numIons - rotRank;
end
else
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if ((tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2) & (rotRank > 1)  
numIons = numIons - 2;
cand = [0 0 tmp(3); 0 0 -1*tmp(3)];             
else      
if abs(tmp(3)*newLattice(3,3)) < minDistance/2
tooClose = 1;
tmp(3) = 0;
else
tooClose = 0;
end 
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
if tooClose == 0
cand1 = cand*Hz;
cand = cat(1,cand,cand1);  
end
if tooClose      
numIons = numIons - rotRank;
else
numIons = numIons - 2*rotRank;
end
end
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if numIons == 0 
break
end
end
elseif  ~isempty(findstr(nsym, 'Dh')) | ~isempty(findstr(nsym, 'dh')) 
rotRank = str2num(nsym(3:end));   
if rotRank > 2
newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
end
candidate = [];
while 1
if numIons == 1;   
tmp = [0 0 0];          
candidate = cat(1,candidate,tmp); 
break;
end
if numIons < 4*rotRank 
if numIons > 2*rotRank      
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm([tmp(1) tmp(2)]) > 0.5)
tmp = rand(1,3) - 0.5;
end
tmp(3) = 0;
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
for i = 1 : rotRank
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
numIons = numIons - 2*rotRank;
elseif numIons < rotRank      
cand = [0 0 rand-0.5];
cand = cat(1,cand,cand*Hz); 
numIons = numIons - 2;
else                   
tmp = rand(1,3) - 0.5; 
r = norm(tmp);
tmp = [r 0 0];  
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
numIons = numIons - rotRank;
end
else
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if ((tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2) & (rotRank > 1)  
numIons = numIons - 2;
cand = [0 0 tmp(3); 0 0 -1*tmp(3)];              
else
if abs(tmp(3)*newLattice(3,3)) < minDistance/2
tooClose1 = 1;
tmp(3) = 0;
else
tooClose1 = 0;
end 
tooClose2 = 0;
for i = 1 : rotRank    
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextAtom = tmp*opemat;
if norm((nextAtom-tmp)*newLattice) < minDistance 
tooClose2 = 1;
tooClose1 = 1;      
tmp = (tmp + nextAtom)/2;
cand = tmp;
break;
end
end
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
if tooClose2 == 0
for i = 1 : rotRank
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
end
if tooClose1 == 0
cand1 = cand*Hz;
cand = cat(1,cand,cand1);  
end
if tooClose2      
numIons = numIons - rotRank;
elseif tooClose1
numIons = numIons - 2*rotRank;
else
numIons = numIons - 4*rotRank;
end
end
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if numIons == 0 
break
end
end
elseif  ~isempty(findstr(nsym, 'Dd')) | ~isempty(findstr(nsym, 'dd')) | ~isempty(findstr(nsym, 'Dv')) | ~isempty(findstr(nsym, 'dv')) 
rotRank = str2num(nsym(3:end));   
if rotRank > 1
newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
end
candidate = [];
while 1
if numIons == 1;   
tmp = [0 0 0];          
candidate = cat(1,candidate,tmp); 
break;
end
if numIons < 4*rotRank 
if numIons < 2*rotRank      
cand = [0 0 rand-0.5];
cand = cat(1,cand,cand*Hz); 
numIons = numIons - 2;
else                   
tmp = rand(1,3) - 0.5; 
r = norm(tmp);
tmp = [r 0 0];  
cand = tmp;
angleM = pi/rotRank;
opematM = [cos(angleM) sin(angleM) 0; sin(angleM) -cos(angleM) 0; 0 0 1];  
mirrorAtom = tmp*opematM;
cand = cat(1,cand,mirrorAtom);
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
nextMirrorAtom = mirrorAtom*opemat;
cand = cat(1,cand,nextMirrorAtom);
end
numIons = numIons - 2*rotRank;
end
else
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  
numIons = numIons - 1;
cand = [0 0 tmp(3)];          
else
cand = tmp;
tooClose1 = 0;
for i = 1 : rotRank        
angle = (2*i*pi/rotRank) + pi/rotRank;    
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
if norm((nextAtom-tmp)*newLattice) < minDistance 
tooClose1 = 1;
tmp = (tmp + nextAtom)/2;
cand = tmp;
break;
end
end
tooClose2 = 0;
for i = 1 : rotRank    
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextAtom = tmp*opemat;
if norm((nextAtom-tmp)*newLattice) < minDistance 
tooClose2 = 1;
tmp = (tmp + nextAtom)/2;
cand = tmp;
break;
end
end
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
if tooClose2 == 0
for i = 1 : rotRank
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
end
if (tooClose2 == 1) | ((tooClose2 == 0) & (tooClose1 == 0))
angleM = pi/rotRank;
opematM = [cos(angleM) sin(angleM) 0; sin(angleM) -cos(angleM) 0; 0 0 1];  
mirrorAtom = tmp*opematM;
cand = cat(1,cand,mirrorAtom);
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextMirrorAtom = mirrorAtom*opemat;
cand = cat(1,cand,nextMirrorAtom);
end
if tooClose2 == 0
for i = 1 : rotRank
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextMirrorAtom = mirrorAtom*opemat;
cand = cat(1,cand,nextMirrorAtom);
end
end
end
if tooClose2 == 1       
numIons = numIons - 2*rotRank;
elseif tooClose1 == 1   
numIons = numIons - 2*rotRank;
else
numIons = numIons - 4*rotRank;
end
end
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if numIons == 0 
break
end
end
elseif  ~isempty(findstr(nsym, 'D')) | ~isempty(findstr(nsym, 'd')) 
rotRank = str2num(nsym(2:end));   
if rotRank > 2
newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
end
candidate = [];
while 1
if numIons == 1;   
tmp = [0 0 0];          
candidate = cat(1,candidate,tmp); 
break;
end
if numIons < 2*rotRank 
if numIons < rotRank      
cand = [0 0 rand-0.5];
cand = cat(1,cand,cand*Hz); 
numIons = numIons - 2;
else                   
tmp = rand(1,3) - 0.5; 
r = norm(tmp);
tmp = [r 0 0];  
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
numIons = numIons - rotRank;
end
else
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  
numIons = numIons - 1;
cand = [0 0 tmp(3)];          
else
cand = tmp;
tooClose = 0;
for i = 1 : rotRank    
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextAtom = tmp*opemat;
if norm((nextAtom-tmp)*newLattice) < minDistance 
tooClose = 1;
tmp = (tmp + nextAtom)/2;
cand = tmp;
break;
end
end
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
if tooClose == 0
for i = 1 : rotRank
angle = 2*i*pi/rotRank;
opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
end
if tooClose
numIons = numIons - rotRank;
else
numIons = numIons - 2*rotRank;
end
end
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if numIons == 0 
break
end
end
elseif ~isempty(findstr(nsym, 'O')) | ~isempty(findstr(nsym, 'o'))   
newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
candidate = [];
if (numIons ~= 1) & (numIons ~= 6) & (numIons ~= 7) & (numIons ~= 8) & (numIons ~= 9) & (numIons < 12)
status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
[nothing, nothing] = unix(['echo ' status ' > error_cluster_symmetry']);
errorS = 1;
end
if mod(numIons,2) == 1 
candidate = [0 0 0];
numIons = numIons - 1;
end
n24 = 0;
numIons1 = numIons;
while 1
numIons1 = numIons1 - 24;
if (numIons1 ~= 6) & (numIons1 ~= 8) & (numIons1 ~= 0) & (numIons1 < 12)
break;
else
n24 = n24 + 1;
end
end
n = 1;
while n <= n24  
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if (strcmp(nsym, 'Oh') | strcmp(nsym, 'oh')) & (mod(n24,2) == 1) & (n == n24)  
w = ceil(rand*2);
if w == 1 
tmp(ceil(rand*3)) = 0;  
else 
w = ceil(rand*3);
if w == 1
tmp(1) = tmp(2);  
elseif w == 2
tmp(1) = tmp(3);  
else
tmp(3) = tmp(2);  
end
end
end
cand = tmp;
for i = 1 : 3
angle = i*pi/2;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
for i = 1 : 3
angle = i*pi/2;
opemat = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
for i = 1 : 3
angle = i*pi/2;
opemat = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
nextAtom = tmp*[0 1 0; 1 0 0; 0 0 -1];
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 -1 0; -1 0 0; 0 0 -1];
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[-1 0 0; 0 0 1; 0 1 0];
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[-1 0 0; 0 0 -1; 0 -1 0];
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 1; 0 -1 0; 1 0 0];
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 -1; 0 -1 0; -1 0 0];
cand = cat(1,cand,nextAtom);    
nextAtom = tmp*[0 1 0; 0 0 1; 1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 1; 1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 1 0; 0 0 -1; -1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 -1; 1 0 0; 0 -1 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 -1 0; 0 0 1; -1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 -1; -1 0 0; 0 1 0]; 
cand = cat(1,cand,nextAtom);     
nextAtom = tmp*[0 -1 0; 0 0 -1; 1 0 0]; 
cand = cat(1,cand,nextAtom);
nextAtom = tmp*[0 0 1; -1 0 0; 0 -1 0]; 
cand = cat(1,cand,nextAtom); 
if (strcmp(nsym, 'Oh') | strcmp(nsym, 'oh')) & ~((mod(n24,2) == 1) & (n == n24))  
numIons = numIons - 48;
n = n + 2;
cand1 = -1*cand;
cand = cat(1,cand,cand1); 
else
numIons = numIons - 24;
n = n + 1;
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if n > n24 
break
end
end
if numIons > 0  
n1 = ceil(numIons/8);
n2 = floor(numIons/6);
n = n1 + round(rand*(n2-n1)); 
x = 4*n - round(numIons/2);
y = round(numIons/2) - 3*n;
for i = 1 : x   
tmp = rand - 0.5;
while abs(tmp) < 0.05
tmp = rand - 0.5;
end
candidate = cat(1, candidate, [0 0 tmp]);
candidate = cat(1, candidate, [0 0 -1*tmp]);
candidate = cat(1, candidate, [0 tmp 0]);
candidate = cat(1, candidate, [0 -1*tmp 0]);
candidate = cat(1, candidate, [tmp 0 0]);
candidate = cat(1, candidate, [-1*tmp 0 0]);
end
for i = 1 : y  
tmp = rand - 0.5;
while abs(tmp) < 0.05 | ((ellipse_mode == 1) & (norm([tmp tmp tmp]) > 0.5))
tmp = rand - 0.5;
end
candidate = cat(1, candidate, [tmp tmp tmp]);
candidate = cat(1, candidate, [tmp tmp -1*tmp]);
candidate = cat(1, candidate, [tmp -1*tmp tmp]);
candidate = cat(1, candidate, [tmp -1*tmp -1*tmp]);
candidate = cat(1, candidate, [-1*tmp tmp tmp]);
candidate = cat(1, candidate, [-1*tmp tmp -1*tmp]);
candidate = cat(1, candidate, [-1*tmp -1*tmp tmp]);
candidate = cat(1, candidate, [-1*tmp -1*tmp -1*tmp]);
end
end
else    
rotRank = str2num(nsym(2:end));   
if rotRank > 2
newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
end
candidate = [];
while 1
if numIons < rotRank    
cand = [0 0 rand-0.5];
for i = 2 : numIons
tmp = [0 0 rand-0.5]; 
cand = cat(1,cand,tmp); 
end
numIons = 0;
else
tmp = rand(1,3) - 0.5; 
while ellipse_mode & (norm(tmp) > 0.5)
tmp = rand(1,3) - 0.5;
end
if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  
numIons = numIons - 1;
cand = [0 0 tmp(3)];          
else
cand = tmp;
for i = 1 : rotRank-1
angle = 2*i*pi/rotRank;
opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  
nextAtom = tmp*opemat;
cand = cat(1,cand,nextAtom);
end
numIons = numIons - rotRank;      
end
end
if isempty(candidate)
candidate = cand; 
else
candidate = cat(1,candidate,cand); 
end
if numIons == 0 
break
end
end
end
AbsoluteCoord = candidate*newLattice;
%unix(['echo cluster > ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo 1.0 >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(1,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(2,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(3,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(nI) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo Direct >> ClusterPOSCAR' num2str(goodPop)]);
%  unix(['echo ' num2str(AbsoluteCoord(i,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
[a,b] = PrincipleAxis(AbsoluteCoord);  
if strcmp(nsym, 'E')
AbsoluteCoord = AbsoluteCoord*a; 
ma = max(AbsoluteCoord);
mi = min(AbsoluteCoord);
newLattice = zeros(3,3);
newLattice(1,1) = ma(1) - mi(1) + 0.02;
newLattice(2,2) = ma(2) - mi(2) + 0.02;
newLattice(3,3) = ma(3) - mi(3) + 0.02;
AbsoluteCoord(:,1) = AbsoluteCoord(:,1) - mi(1) + 0.01; 
AbsoluteCoord(:,2) = AbsoluteCoord(:,2) - mi(2) + 0.01;
AbsoluteCoord(:,3) = AbsoluteCoord(:,3) - mi(3) + 0.01;
candidate = (AbsoluteCoord/newLattice) - 0.5;
elseif ~isempty(findstr(nsym, 'O')) | ~isempty(findstr(nsym, 'o')) | ~isempty(findstr(nsym, 'T')) | ~isempty(findstr(nsym, 't')) 
candidate = candidate*a; 
m = max(max(abs(candidate)));
candidate = candidate*(0.5/(m + 0.01));
else 
AbsoluteCoord = AbsoluteCoord*a; 
mass_center = zeros(1,3);
for i = 1 : 3 
mass_center(i) = sum(AbsoluteCoord(:,i))/nI;
AbsoluteCoord(:,i) = AbsoluteCoord(:,i) - mass_center(i); 
end
m = max(abs(AbsoluteCoord));
newLattice = zeros(3,3);
newLattice(1,1) = 2 * m(1) + 0.02;
newLattice(2,2) = 2 * m(2) + 0.02;
newLattice(3,3) = 2 * m(3) + 0.02;
candidate = AbsoluteCoord/newLattice;
end
candidate = candidate + 0.5;  
%unix(['echo cluster >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo 1.0 >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(1,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(2,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(3,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(nI) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo Direct >> ClusterPOSCAR' num2str(goodPop)]);
%  unix(['echo ' num2str(candidate(i,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
