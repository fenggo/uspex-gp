function [candidate, Lattice, errorS] = symope_crystal(nsym, numIons, lat, minD, sym_coef)
global ORG_STRUC
fixLat = ORG_STRUC.constLattice;    
spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];   
fixRndSeed = ORG_STRUC.fixRndSeed;                              
sGroup = spaceGroups(nsym);                                     
candidate = zeros(sum(numIons),3);                              
Lattice = [1 0 0; 0 1 0; 0 0 1];                                
errorS = 0;                                                     
[permutation, permutationBack] = GetPermutation(nsym);                        
Init_lat                       = Get_Init_Lattice(lat, permutation);          
[Init_numIons, Init_lat,fail] = GetPrimitiveCell(sGroup, Init_lat, numIons,fixLat);
if ~fail
Write_Stokes_input(Init_numIons, minD, nsym, Init_lat, fixRndSeed, sym_coef); 
command = [spgBINDIR '/random_cell < rc.in > rc.out'];                        
[a, b] = unix(command);
[coordinate_S, lattice_S, failed]=Read_Stokes_output('rc.out', Init_numIons); 
if ~failed                                                                    
[lattice, coordinate] = fix_latticeStokes_after(nsym, lattice_S, coordinate_S);
lattice(4:6) = lattice(4:6)*pi/180;
if fixLat == 1
Lattice_Matrix = latConverter(lattice);
else
Lattice_Matrix = latConverter(lattice);
abs_cand = coordinate*Lattice_Matrix;
[abs_cand,Lattice_Matrix] = optLattice(abs_cand, Lattice_Matrix);      
coordinate = abs_cand/Lattice_Matrix;
end
Lattice_Matrix(find(abs(Lattice_Matrix)<0.000001))=0;
coordinate(find(abs(coordinate  )<0.000001))=0;
coordinate(find(abs(coordinate-1)<0.000001))=0;
if latticeCheck(Lattice_Matrix)                                           
Lat_type = sGroup(1);
if (Lat_type ~= 'P') & (fixLat == 1)
if (nsym > 37) & (nsym < 42)                                      
Lat_type = 'B';
end
[Lattice_Matrix, candidate] = unitCellFromPrimitive(Lat_type, Lattice_Matrix, coordinate);
else
candidate = coordinate;
end                                                                   
[Lattice, candidate] = Get_Final_Struc(Lattice_Matrix, candidate, permutationBack);
else
errorS = 1;
end
if size(candidate, 1) ~= sum(numIons)
errorS = 1;
end
else
errorS = 1;
end
end
