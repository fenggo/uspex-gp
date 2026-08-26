function USPEXmessage(msgID, msgString, option)
global ORG_STRUC
clockTIME=clock;
timeString = ['===== USPEX WARNING @ ',num2str(clockTIME(4)),':',num2str(clockTIME(5)),':',num2str(clockTIME(3)), ' ', date];
if msgID > 0	
displayMsg = [ timeString, ' ===== \n', messageInfo(msgID)];
else
displayMsg = [ timeString, ' ===== \n'];
end
if length(msgString)>0
displayMsg = [displayMsg, msgString, '\n\n'];
else
displayMsg = [displayMsg, '\n\n'];
end
if option <= 0
fprintf(1, displayMsg);
end
if option >=0
filePath=[ORG_STRUC.homePath, '/Warnings'];
fp = fopen(filePath, 'a+');
fprintf(fp, displayMsg);
fclose(fp);
end
