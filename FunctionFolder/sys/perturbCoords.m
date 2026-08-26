function Coords = perturbCoords(Coordinates, lattice, pertRegime)
pert_rate = 0.05; 
N = size(Coordinates, 1);
perturbation = 2*(rand(N,3)-0.5)*pert_rate;
if pertRegime == 2   
perturbation(:,1) = perturbation(1,1);
perturbation(:,2) = perturbation(1,2);
perturbation(:,3) = perturbation(1,3);
elseif pertRegime == 3    
perturbation(:,3) = 0;
end
perturbation = perturbation/lattice; 
Coords = Coordinates + perturbation;
if pertRegime ~= 2
Coords = Coords - floor(Coords);
end