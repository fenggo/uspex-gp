function  dx_pro = NEB_trustRadius(dX, N)
global ORG_STRUC
global POP_STRUC
trustRadius = 0.25;     
maxRelaxVolume = 10/100;   
sumIons  = sum(ORG_STRUC.numIons);
dimension  = ORG_STRUC.dimension;
numDimension = 3*(sumIons+dimension);
dx_pro = zeros(numDimension,1);
dx_pro=dX;
tmp = ( reshape( dX',3,sumIons+dimension ) )';
CARTECOORDSDisplace = tmp(1+dimension:sumIons+dimension, :);
LATTICEStrain    = tmp(1:dimension,:);
epsil = LATTICEStrain;
if trace(epsil) > maxRelaxVolume
dx_pro(1:9) = dX(1:9) * maxRelaxVolume/trace(epsil);
%   disp('displacement of the lattic cell :')
end
displacement = zeros(1, sumIons);
for i = 1:sumIons
displacement(i) = norm( CARTECOORDSDisplace(i,:) ) ;
end
if max( displacement ) > trustRadius
dx_pro(10:end) = dX(10:end) * trustRadius/max( displacement ); 
end
