function density = calcDensity( numIons, atomType, Volume)
atomicMass=zeros(1, length(atomType));
for i = 1:length(atomType)
atomicMass(i) = elementMass(atomType(i));
end
density = 'N/A';
if ~strcmp(Volume, 'N/A')
density = sum(numIons.*atomicMass)/Volume/0.602214129;
end
