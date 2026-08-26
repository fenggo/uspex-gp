function lat_OK = latticeCheck(Lattice)
global ORG_STRUC
lat_OK = 1;
if length(Lattice)==3
angLattice = latConverter(Lattice);
else
angLattice = Lattice;
Lattice = latConverter(Lattice);
end
angles = angLattice(4:6)*180/pi;
if (isempty(find(angles<ORG_STRUC.minAngle)) + isempty(find(angles>(180-ORG_STRUC.minAngle))))~=2
lat_OK = 0;
end
vol = abs(det(Lattice));
dist = zeros(3,1);
dist(1) = vol/(angLattice(1)*angLattice(2)*sin(angLattice(6)));
dist(2) = vol/(angLattice(1)*angLattice(3)*sin(angLattice(5)));
dist(3) = vol/(angLattice(2)*angLattice(3)*sin(angLattice(4)));
if ~isempty(find(angLattice(1:3)<ORG_STRUC.minVectorLength))
lat_OK = 0;
end
if ~isreal(Lattice) || sum(find(isnan(Lattice))) > 0  
lat_OK = 0;
end
if lat_OK 
a_bc = Lattice(1,1)*(Lattice(2,1)+Lattice(3,1))+Lattice(1,2)*(Lattice(2,2)+Lattice(3,2))+Lattice(1,3)*(Lattice(2,3)+Lattice(3,3));
ab_c = Lattice(3,1)*(Lattice(2,1)+Lattice(1,1))+Lattice(3,2)*(Lattice(2,2)+Lattice(1,2))+Lattice(3,3)*(Lattice(2,3)+Lattice(1,3));
b_ca = Lattice(2,1)*(Lattice(1,1)+Lattice(3,1))+Lattice(2,2)*(Lattice(1,2)+Lattice(3,2))+Lattice(2,3)*(Lattice(1,3)+Lattice(3,3));
anglesDiag1(1) = acos(a_bc/(angLattice(1)*sqrt(angLattice(2)^2+angLattice(3)^2+2*angLattice(2)*angLattice(3)*cos(angLattice(4)))));
anglesDiag1(2) = acos(b_ca/(angLattice(2)*sqrt(angLattice(1)^2+angLattice(3)^2+2*angLattice(1)*angLattice(3)*cos(angLattice(5)))));
anglesDiag1(3) = acos(ab_c/(angLattice(3)*sqrt(angLattice(2)^2+angLattice(1)^2+2*angLattice(2)*angLattice(1)*cos(angLattice(6)))));
anglesDiag = anglesDiag1(1:3)*180/pi;
if (isempty(find(anglesDiag<ORG_STRUC.minDiagAngle)) + isempty(find(anglesDiag>(180-ORG_STRUC.minDiagAngle))))~=2
lat_OK = 0;
end
end
if ORG_STRUC.dimension==2
lat_OK = 1
end
if ORG_STRUC.constLattice
lat_OK = 1;
end
