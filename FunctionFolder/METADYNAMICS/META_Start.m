function META_Start()
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
rand('state',sum(100*clock));
c = clock;
rand_init = round(c(5)+c(6));
for i = 1:rand_init
dummy = rand;
dummy = randperm(5);
end
disp('  ');
while 1
stillRunningInfo();
[nothing,homePath] = unix('pwd');
homePath(end) = [];
USPEXPath= homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
if findstr(uspexmode,'exe')
[nothing,USPEXPath]= unix('echo -e $USPEXPATH');;
USPEXPath(end)=[];
end
if exist('Current_POP.mat')
load Current_POP.mat
load Current_ORG.mat
if ~(findstr(ORG_STRUC.homePath, homePath)) || length(ORG_STRUC.homePath)~=length(homePath)
mesg=[ORG_STRUC.homePath ' --> ' homePath ];
ORG_STRUC.homePath=homePath;
USPEXmessage(1011, mesg, 0);
end
ORG_STRUC.USPEXPath=USPEXPath;
cd (POP_STRUC.resFolder)
load USPEX.mat
cd ..
else
inputFile = 'INPUT.txt';
META_CreateORGStruc(inputFile);
[nothing, nothing] = unix (['cp ' inputFile ' ' ORG_STRUC.resFolder '/Parameters.txt']);
ORG_STRUC.pickedUP = 0;
CreateCalcFolder(1);
if ORG_STRUC.pickUpYN
PickUp();
else
META_Initialize();
end
CreateCalcFolder();
[nothing, nothing] = unix (['cp POSCAR_1 ' POP_STRUC.resFolder '/']);
end
safesave ('Current_POP.mat', POP_STRUC)
safesave ('Current_ORG.mat', ORG_STRUC)
META();
if ORG_STRUC.repeatForStatistics > 1
RepeatRun();
else
if exist('still_running')
[nothing, nothing] = unix('rm still_running');
[nothing, nothing] = unix('rm POSCAR');
[nothing, nothing] = unix('rm POSCAR_order');
end
disp('USPEX IS DONE!');
[nothing, nothing] = unix('echo 1 > USPEX_IS_DONE');
disp(' ');
quit
end
end
