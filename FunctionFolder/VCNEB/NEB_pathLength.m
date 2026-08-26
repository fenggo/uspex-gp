function NEB_pathLength()
global ORG_STRUC
global POP_STRUC
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);
dimension  = ORG_STRUC.dimension;
numDimension = 3*(sumIons+dimension);
for i=1:numImages
L = POP_STRUC.POPULATION(i).LATTICE;
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-2)-2) = L(1,1);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-2)-1) = L(1,2);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-2)-0) = L(1,3);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-1)-2) = L(2,1);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-1)-1) = L(2,2);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-1)-0) = L(2,3);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-0)-2) = L(3,1);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-0)-1) = L(3,2);
POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-0)-0) = L(3,3);
U = POP_STRUC.POPULATION(i).caCoordVector(1:3*dimension);
POP_STRUC.POPULATION(i).frCoordVector(1:3) = U(1:3)/norm( U(1:3));
POP_STRUC.POPULATION(i).frCoordVector(4:6) = U(4:6)/norm( U(4:6));
POP_STRUC.POPULATION(i).frCoordVector(7:9) = U(7:9)/norm( U(7:9));
for j = 1+dimension:sumIons+dimension    
sx = POP_STRUC.POPULATION(i).CARTECOORDS(j-dimension,1);   
sy = POP_STRUC.POPULATION(i).CARTECOORDS(j-dimension,2);   
sz = POP_STRUC.POPULATION(i).CARTECOORDS(j-dimension,3);   
POP_STRUC.POPULATION(i).caCoordVector( (j-1)*3+1:j*3) = [sx, sy, sz];
end 
for j = 1+dimension:sumIons+dimension    
sx = POP_STRUC.POPULATION(i).COORDINATES(j-dimension,1);   
sy = POP_STRUC.POPULATION(i).COORDINATES(j-dimension,2);   
sz = POP_STRUC.POPULATION(i).COORDINATES(j-dimension,3);   
POP_STRUC.POPULATION(i).frCoordVector( (j-1)*3+1:j*3) = [sx, sy, sz];
end
end
pathLength=zeros(1, numImages);
%disp('PATH length')
POP_STRUC.POPULATION(1).pathLength = 0;
for i = 2:numImages
POP_STRUC.POPULATION(i).pathLength  =  norm( POP_STRUC.POPULATION(i).caCoordVector - POP_STRUC.POPULATION(i-1).caCoordVector );
pathLength(i) = POP_STRUC.POPULATION(i).pathLength ;
end
