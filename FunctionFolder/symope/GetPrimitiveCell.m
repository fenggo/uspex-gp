function [numIons, lat1, Error] = GetPrimitiveCell(sGroup, lat1, numIons, fixLat);
Error = 0;
if (sGroup(1) ~= 'P') & (fixLat == 1)            
numIons0 = numIons;
lat0     = lat1;
if sGroup(1) == 'F'                    
Duplicate = 4;      
elseif sGroup(1) == 'R'
[doedl, IX] = sort(lat1(4:6));
if ( norm(doedl'-[90 90 120]) < 3 )       
Duplicate = 3;   
else
Duplicate = 1;                       
lat2 = latConverter(lat1);           
lat2 = [1 -1 0; 0 1 -1; 1 1 1]*lat2; 
lat1 = latConverter(lat2);
lat1(4:6) = lat1(4:6)*180/pi; 
end
else
Duplicate = 2;                           
end
numIons = round(numIons/Duplicate);
lat1(1:3) = lat1(1:3)/power(Duplicate,1/3);
if mod(sum(numIons0),sum(numIons))>0         
disp(['impossible to place ' num2str(numIons) ' for ' sGroup ' symmetry']);
Error = 1;
end
end
