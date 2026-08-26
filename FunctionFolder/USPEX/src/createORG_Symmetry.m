function createORG_Symmetry(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
doSpaceGroup = python_uspex(getPy, ['-f ' inputFile ' -b doSpaceGroup -c 1']);
if isempty(doSpaceGroup)
if (ORG_STRUC.dimension==0) | (ORG_STRUC.dimension==2) | (ORG_STRUC.dimension==-3)
doSpaceGroup = '0'; 
else
doSpaceGroup = '1'; 
end
end
ORG_STRUC.doSpaceGroup = str2num(doSpaceGroup);
symmetrize = python_uspex(getPy, ['-f ' inputFile ' -b symmetrize -c 1']);
if ~isempty(symmetrize)
ORG_STRUC.symmetrize = str2num(symmetrize);
end
SGtolerance = python_uspex(getPy, ['-f ' inputFile ' -b SymTolerance -c 1']);
if ~isempty(SGtolerance)
SGtolerance = deblank(SGtolerance);
if strcmp(lower(SGtolerance), 'high')
ORG_STRUC.SGtolerance = 0.05;
elseif strcmp(lower(SGtolerance), 'medium')
ORG_STRUC.SGtolerance = 0.10;
elseif strcmp(lower(SGtolerance), 'low')
ORG_STRUC.SGtolerance = 0.20;
else 
ORG_STRUC.SGtolerance = str2num(SGtolerance);
end
end
if isempty(ORG_STRUC.SGtolerance) 
ORG_STRUC.SGtolerance = 0.10;
end
nsym = python_uspex(getPy, ['-f ' inputFile ' -b symmetries -e EndSymmetries']);
if isempty(nsym)
if ORG_STRUC.dimension==0
nsym = 'E C2 D2 C4 C3 C6 T S2 Ch1 Cv2 S4 S6 Ch3 Th Ch2 Dh2 Ch4 D3 Ch6 O D4 Cv3 D6 Td Cv4 Dd3 Cv6 Oh Dd2 Dh3 Dh4 Dh6 Oh C5 S5 S10 Cv5 Ch5 D5 Dd5 Dh5 I Ih '; 
elseif ORG_STRUC.dimension==3
nsym = '2-230'; 
elseif ORG_STRUC.dimension==-2 | ORG_STRUC.dimension==1 | ORG_STRUC.dimension==-3
nsym = '2-17';
end
end
if ORG_STRUC.dimension==0
nsym(end) = [];
ORG_STRUC.nsym = nsym;
c1 = findstr(nsym, ' ');
c = sort(str2num(['0 ' num2str(c1)]));
c(end+1) = length(nsym) + 1;
ind = zeros(1,2);
indN = 0;
for i = 2 : length(c)
if c(i-1)+1 > c(i)-1
continue
end
indN = indN + 1;
ind(1,1) = c(i-1)+1;
ind(1,2) = c(i)-1;
if indN == 1
ORG_STRUC.nsymN = ind;
else
ORG_STRUC.nsymN = cat(1,ORG_STRUC.nsymN,ind);
end
end
if isempty(ORG_STRUC.nsym)
ORG_STRUC.nsym = 'E';
ORG_STRUC.nsymN = [1 1];
end
else
ORG_STRUC.nsym = zeros(1,230);
c1 = findstr(nsym, ' ');
c2 = findstr(nsym, '-');
c = sort(str2num(['0 ' num2str(c1) ' ' num2str(c2)]));
c(end+1) = length(nsym) + 1;
ind1 = 1;
for i = 2 : length(c)
if c(i-1)+1 > c(i)-1
continue
end
ind2 = str2num(nsym(c(i-1)+1 : c(i)-1));
if ind2 == 0
ind1 = 1;
continue
end
if ~isempty(find(c2 == c(i-1)))
for j = ind1 : ind2
ORG_STRUC.nsym(j) = 1;
end
else
ORG_STRUC.nsym(ind2) = 1;
end
ind1 = ind2;
end
if sum(ORG_STRUC.nsym) == 0
ORG_STRUC.nsym(1) = 1;
end
ORG_STRUC.nsymN = [0 0];
end
sym_coef = python_uspex(getPy, ['-f ' inputFile ' -b constraint_enhancement -c 1']);
if ~isempty(sym_coef)
ORG_STRUC.sym_coef = str2num(sym_coef);
end
