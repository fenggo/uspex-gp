function [lat] = fix_latticeStokes_before(nsym, oldLat)
lat = oldLat;
sGroup = spaceGroups(nsym); 
if (nsym > 2) & (nsym < 16) & (sGroup(1) == 'P') 
lat(1) = oldLat(3); 
lat(3) = oldLat(1); 
end
