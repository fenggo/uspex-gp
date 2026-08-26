function [lat, coord] = fix_latticeStokes_after(nsym, oldLat, oldCoord)
lat = oldLat;
coord = oldCoord;
sGroup = spaceGroups(nsym); 
if (nsym > 2) & (nsym < 16) & (sGroup(1) == 'P') 
lat(5) = oldLat(6); 
lat(6) = 90;
lat(4) = 90;
lat(1) = oldLat(2); 
lat(2) = oldLat(3); 
lat(3) = oldLat(1); 
coord(:,1) = oldCoord(:,2); 
coord(:,2) = oldCoord(:,3);
coord(:,3) = oldCoord(:,1);
end
