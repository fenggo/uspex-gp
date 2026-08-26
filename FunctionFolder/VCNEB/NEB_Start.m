function NEB_Start()
global ORG_STRUC
global POP_STRUC
disp('  ');
while 1
if exist('DONE_VCNEB') 
quit();
end
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
if isempty(findstr(ORG_STRUC.homePath, homePath)) || length(ORG_STRUC.homePath)~=length(homePath)
mesg=[ORG_STRUC.homePath ' --> ' homePath ];
ORG_STRUC.homePath=homePath;
USPEXmessage(1011, mesg, 0);
end
ORG_STRUC.USPEXPath=USPEXPath;
cd (POP_STRUC.resFolder)
else
if exist('INPUT.txt')
inputFile = 'INPUT.txt';
else
inputFile = 'INPUT_VCNEB.txt';
end
createORGStruc(inputFile);
[nothing, nothing] = unix (['cp ' inputFile ' ' ORG_STRUC.resFolder '/Parameters.txt']);
[nothing, nothing] = unix (['cp Images ' ORG_STRUC.resFolder '/initialImages']);
ORG_STRUC.pickedUP = 0;
CreateCalcFolder(1);
if ORG_STRUC.pickUpYN
NEB_PickUp();
else
NEB_initialize_POP_STRUC();
end
CreateCalcFolder();
POP_STRUC.resFolder = ORG_STRUC.resFolder;
end
safesave ([ORG_STRUC.homePath, '/Current_POP.mat'], POP_STRUC)
safesave ([ORG_STRUC.homePath, '/Current_ORG.mat'], ORG_STRUC)
VCNEB();
if ORG_STRUC.repeatForStatistics > 1
RepeatRun();
else
if exist('still_running')
[nothing, nothing] = unix('rm still_running');
end
disp('VCNEB IS DONE!');
[nothing, nothing] = unix('echo 1 > VCNEB_IS_DONE');
disp(' ');
quit
end
end
