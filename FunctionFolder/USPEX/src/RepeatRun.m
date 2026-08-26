function RepeatRun()
global ORG_STRUC
flag = 1;
rfs=python_uspex([ORG_STRUC.USPEXPath '/FunctionFolder/getInput.py'], ['-f multiple_runs -b runs_left -c 1']);
rfs_1 = str2num(rfs);
if rfs_1 == 1
disp('USPEX IS DONE!');
[nothing, nothing] = unix('echo 1 > USPEX_IS_DONE');
disp(' ');
quit
end
if (rfs_1 > 1) && (flag)
rfs_1 = rfs_1 - 1;
rfs = [num2str(rfs_1) ' : runs_left'];
[nothing, nothing] = unix(['echo "' rfs '" > multiple_runs']);
[nothing, nothing] = unix('rm -rf CalcFold*');
[nothing, nothing] = unix('rm -rf *.mat');
end
