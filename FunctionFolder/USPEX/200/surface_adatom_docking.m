function [ z_fractional ] = surface_adatom_docking( lattice, zmin, substrate_coords, substrate_atom_type, adatom_coord, adatom_type )
z_fractional = zmin ;
z_cart = z_fractional * lattice(3,3) ;
natoms = numel( substrate_coords(:,1) );
for i = 1:natoms 
sub_z = lattice(3,3) * substrate_coords(i,3) ;
parallel_dist = distance( lattice, substrate_coords(i,:), [ adatom_coord, substrate_coords(i,3) ] );
radii1 = str2num(covalentRadius(ceil(substrate_atom_type(i))));
radii2 = str2num(covalentRadius(ceil(adatom_type)));
min_dist = 1.0*(radii1+radii2);
if parallel_dist < min_dist
z_shift = sqrt( min_dist^2 - parallel_dist^2 ) ;
if z_cart < z_shift + sub_z 
z_cart = z_shift + sub_z;
z_fractional = z_cart / lattice(3,3) ;
end
end
end
function dist=distance(lattice, coord1, coord2)
check = coord1-coord2;
for index=1:3
check(index)=check(index)-round(check(index));
end
dist = sqrt(sum((check*lattice).^2));
