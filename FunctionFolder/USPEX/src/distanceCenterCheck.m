function constraintsOK = distanceCenterCheck (Coordinates,Lattice,composition)
global ORG_STRUC
dist =[];
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
minimalDistance=10;
for loop = 1:(size(Coordinates,1)-1)
if constraintsOK    
for who = (loop+1):size(Coordinates,1)
row = find(ionChange>(loop-0.5));
colomn = find(ionChange>(who-0.5));
check = (Coordinates(loop,:)-Coordinates(who,:));
for x=1:3
for y=1:3
for z=1:3
dist = sqrt(sum(((check+[x-2,y-2,z-2])*Lattice).^2));
if dist<minimalDistance
minimalDistance = dist;
end
if dist < (ORG_STRUC.CenterminDistMatrice(row(1),colomn(1))-0.1)
constraintsOK = 0;
break   
end
end
end
end
row =[];
colomn = [];
dist = [];
check = [];
where = [];
end
else                
break           
end                 
end
