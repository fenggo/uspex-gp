function constraintsOK = distanceCheck(Coordinates, Lattice, composition, minDistMatrice)
global ORG_STRUC
if isempty(find(minDistMatrice>0.0001))
constraintsOK = 1
elseif size(Coordinates,1) ~= sum(composition)
constraintsOK = 0
disp(['Opps: Size of coordinates <> sum(numIons)!']);
debug1 = Coordinates
debug2 = composition
elseif ~isreal(Coordinates) || ~isreal(Lattice)
constraintsOK = 0
disp(['Opps: complex number appears in coordinates/lattice!']);
debug1 = Coordinates
debug2 = Lattice
elseif sum(find(isnan(Coordinates))) > 0 || sum(find(isnan(Lattice))) > 0
constraintsOK = 0;
% disp(['Opps: NaN appears in coordinates/lattice!']);
debug1 = Lattice;
debug2 = Coordinates;
else
dist = [];
constraintsOK = 1;
breakCounter = 0;
ionChange = zeros(1,length(composition));
for ind = 1:length(composition)
ionChange(ind) = sum(composition(1:ind));
end
if length(Lattice)==6
Lattice = latConverter(Lattice);
end
if ~isempty(find(size(Coordinates)==1))
Coordinates_temp = Coordinates;
Coordinates = zeros(3,sum(composition));
Coordinates(:) = Coordinates_temp(:);
Coordinates = Coordinates';
end
if ORG_STRUC.dimension == 0
Matrix = [0 0 0];
else
X1 = [-1:1];
Y1 = [-1:1];
Z1 = [-1:1];
[X2,Y2,Z2] = meshgrid(X1,Y1,Z1);
X3=reshape(X2,1,[]);
Y3=reshape(Y2,1,[]);
Z3=reshape(Z2,1,[]);
Matrix = [X3; Y3; Z3]';
end
for loop = 1:(size(Coordinates,1)-1)
if constraintsOK    
for who = (loop+1):size(Coordinates,1)
row = find(ionChange>(loop-0.5));
colomn = find(ionChange>(who-0.5));
coor1 = Coordinates(loop,:)*Lattice;
coor2 = Coordinates(who,:)*Lattice;
tmp   = bsxfun(@plus, Matrix*Lattice, coor2);
dist  = min(pdist2(tmp, coor1));
if dist < minDistMatrice(row(1),colomn(1))
constraintsOK = 0;
break   
end
end
else                
break           
end                 
end
end
