function NEB_initialize_POP_STRUC()
global ORG_STRUC
global POP_STRUC
POP_STRUC = NEB_creatPOPStruct(); 
numImages = ORG_STRUC.numImages;
numIons = ORG_STRUC.numIons;
dimension = ORG_STRUC.dimension;
numDimension = 3*( sum(numIons) + dimension );
POP_STRUC.genDone = 0;
POP_STRUC.step = 0;
POP_STRUC.bodyCount = 0;
POP_STRUC.DoneOrder = zeros(numImages,1);
if ~exist('Images')
disp('  ');
disp('VCNEB ERROR : File "Images" for VCNEB is not found ! ...');
disp('  ');
exit
end
fpath = [ORG_STRUC.homePath,'/Images'];
fp = fopen(fpath);
[nothing, maxImageINPUT] = unix(['grep Image ' fpath ' | wc -l' ]);
for i = 1:str2num(maxImageINPUT)
tline = fgetl(fp);
if      (findstr(tline, 'Image_ini'))
pos = 1;
elseif  (findstr(tline, 'Image_end'))
pos = numImages;
elseif  (findstr(tline, 'Image '))
pos = sscanf(tline, 'Image %d');
else
break;
end
scale_factor = str2num(fgetl(fp));
lat = fscanf(fp, '%g\n', [3,3]);
if pos > numImages
else
POP_STRUC.POPULATION(pos).LATTICE = lat'*scale_factor;
end
num = str2num(fgetl(fp));
if sum(num) ~= sum(numIons)
disp('VCNEB ERROR : The number of atoms is inconsistent from INPUT.txt and Image');
disp('VCNEB ERROR : VCNEB JOB QUIT WITH ERROR !')
quit
end
tmp = fgetl(fp);
if pos > numImages
for j = 1:sum(num)
str2num(fgetl(fp));
end
else
for j = 1:sum(num)
POP_STRUC.POPULATION(pos).COORDINATES(j,:) = str2num(fgetl(fp));
end
POP_STRUC.POPULATION(pos).CARTECOORDS = POP_STRUC.POPULATION(pos).COORDINATES * POP_STRUC.POPULATION(pos).LATTICE;
end
end
%disp(Atom.ForSymfind);
Vini= abs( det( POP_STRUC.POPULATION(1).LATTICE ) );
Vfin= abs( det( POP_STRUC.POPULATION(numImages).LATTICE ) );
if 0==Vini
disp('VCNEB ERR: Wrong setting in 1-st Image! Please check the Image file.');
quit
end
if 0==Vfin
disp('VCNEB ERR: Wrong setting in the Final Image! Please check the Image file.');
quit
end
if       0==ORG_STRUC.optReadImages
disp('VCNEB MSG: Reading all the images from the Image file !');
elseif   1==ORG_STRUC.optReadImages
disp('VCNEB MSG: Using the initial and filanl Images to generate the Images !');
for i = 2:numImages-1
tmp1 = zeros(sum(numIons),3);
tmp2 = zeros(3,3);
tmp1 = (1-i)/(numImages-1)*(POP_STRUC.POPULATION(1).CARTECOORDS - POP_STRUC.POPULATION(numImages).CARTECOORDS);
POP_STRUC.POPULATION(i).CARTECOORDS = tmp1 + POP_STRUC.POPULATION(1).CARTECOORDS;
tmp2 = (1-i)/(numImages-1)*(POP_STRUC.POPULATION(1).LATTICE - POP_STRUC.POPULATION(numImages).LATTICE);
POP_STRUC.POPULATION(i).LATTICE = tmp2 + POP_STRUC.POPULATION(1).LATTICE;
end
elseif  2==ORG_STRUC.optReadImages
disp('VCNEB MSG: Automatically generate the Images based on the Image file!')
allImageVolume = zeros(numImages,1);
for i = 1:numImages
allImageVolume(i)= det( POP_STRUC.POPULATION(i).LATTICE );
end
for i = 2:numImages
if allImageVolume(i) == 0
prevImage=i-1; nextImage=0;
for j=i+1:1:numImages
if allImageVolume(j)>0
nextImage=j;
break;
end
end 
diffN = nextImage - prevImage;
D_CARTECOORDS=POP_STRUC.POPULATION(nextImage).CARTECOORDS-POP_STRUC.POPULATION(prevImage).CARTECOORDS;
D.LATTICE   =POP_STRUC.POPULATION(nextImage).LATTICE   -POP_STRUC.POPULATION(prevImage).LATTICE;
POP_STRUC.POPULATION(i).CARTECOORDS=POP_STRUC.POPULATION(prevImage).CARTECOORDS+ D_CARTECOORDS/diffN; 
POP_STRUC.POPULATION(i).LATTICE   =POP_STRUC.POPULATION(prevImage).LATTICE   +    D.LATTICE/diffN;
end
end
end
disp('VCNEB MSG: Checking the Images ...')
for i = 2:numImages-1
if  0 == abs( det( POP_STRUC.POPULATION(i).LATTICE ) )
disp(['VCNEB ERR: The Image' num2str(i) ' error with a VOLUME=0, please check the INPUT.txt and Image files !'])
disp('VCNEB ERR: VCNEB JOB QUIT WITH ERROR !')
quit;
end
end
for i=1:numImages
POP_STRUC.POPULATION(i).COORDINATES = (  ( (POP_STRUC.POPULATION(i).LATTICE')^(-1) )*POP_STRUC.POPULATION(i).CARTECOORDS' )';
if 2 == ORG_STRUC.UnitType
POP_STRUC.POPULATION(i).LATTICE = POP_STRUC.POPULATION(i).LATTICE * ORG_STRUC.bohr2A;
end
end
if     1==ORG_STRUC.FormatType   
XSF_WriteOutFile;
elseif 2==ORG_STRUC.FormatType | 0==ORG_STRUC.FormatType
STM4_WriteOutFile(0);
elseif 3==ORG_STRUC.FormatType
XYZ_WriteOutFile; 
end
disp('VCNEB MSG: Images are successfully generated !')
disp('VCNEB MSG: Start the main code ... ')
disp('');
NEB_pathLength();
if ORG_STRUC.VarPathLength==0 & ORG_STRUC.CalcType==1
for i = 2:numImages
ORG_STRUC.VarPathLength=ORG_STRUC.VarPathLength+POP_STRUC.POPULATION(i).pathLength;
end
ORG_STRUC.VarPathLength=ORG_STRUC.VarPathLength/(numImages-1);
end
NEB_writeOutput(-1);
