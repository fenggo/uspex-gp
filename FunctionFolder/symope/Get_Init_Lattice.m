function lat1 = Get_Init_Lattice(lat, permutation)
if (size(lat) == [1 1])        
lat1 = zeros(1,6);
lat1(1:3) = rand(1,3) + 0.5;
x = 0;
while (x < 0.3)
lat1(4:6) = (rand(1,3)*120 + 30)*pi/180;
x = 1 - cos(lat1(4))^2 - cos(lat1(5))^2 - cos(lat1(6))^2 + 2*cos(lat1(4))*cos(lat1(5))*cos(lat1(6));
end
ratio = lat/det(latConverter(lat1));
lat1(4:6) = lat1(4:6)*180/pi; 
lat1(1:3) = lat1(1:3)*(ratio)^(1/3);
else                           
if (size(lat,1) == 3) & (size(lat,2) == 3) 
lat1 = latConverter(lat);
else                                       
lat1 = lat;
end
lat_tmp = lat1;
for axis = 1 : 3
lat1(axis) = lat_tmp(permutation(axis));
end
lat1(4:6) = lat1(4:6)*180/pi; 
end
