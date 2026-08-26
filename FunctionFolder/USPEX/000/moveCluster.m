function [candidate] = moveCluster(lat,coordinates)
for j = 1 : 3
badPositions = find(coordinates(:,j) > 0.9);
i=0;
while ~isempty(badPositions)
coordinates(:,j) = coordinates(:,j) + (1-coordinates(badPositions(1),j)) + 0.01;
coordinates = coordinates - floor(coordinates);
badPositions = find(coordinates(:,j) > 0.9);
i = i + 1;
if i > 10
%disp('we have a problem here, cluster atoms are all over the cell');
break   
end
end
end
AbsoluteCoord = coordinates*lat;
AbsoluteCoord = bsxfun(@minus, AbsoluteCoord, mean(AbsoluteCoord)); 
[a,b] = PrincipleAxis(AbsoluteCoord);  
candidate = bsxfun(@plus, AbsoluteCoord*a/lat, [0.5, 0.5, 0.5]); 
