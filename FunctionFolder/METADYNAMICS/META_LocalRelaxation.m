function META_LocalRelaxation()
global ORG_STRUC
global POP_STRUC
while 1 
META_ReadJobs();
META_SubmitJobs();
if ORG_STRUC.platform > 0
if sum([POP_STRUC.POPULATION(:).Done]) ~= length(POP_STRUC.POPULATION)
[nothing, nothing] = unix('rm still_running');
fclose ('all');
quit
else
break; 
end
else
if sum([POP_STRUC.POPULATION(:).Done]) == length(POP_STRUC.POPULATION)
break; 
end
end
end   
