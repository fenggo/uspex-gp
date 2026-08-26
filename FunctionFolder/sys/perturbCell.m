function lat = perturbCell(lattice, pertRegime)
global ORG_STRUC
pert_rate_cell = 0.5;
lat = lattice; 
if size(ORG_STRUC.lattice,1) ~= 3   
N = size(lattice, 1);
perturbation = (rand(N,3)-0.5)*pert_rate_cell;
lat = lattice + perturbation;
end
