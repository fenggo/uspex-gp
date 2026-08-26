function RepeatRun()
flag = 1;
[nothing, rfs] = unix (['./getStuff multiple_runs runs_left 1']);
rfs_1 = str2num(rfs);
if rfs_1 == 1
quit
end
if (rfs_1 > 1) && (flag)
rfs_1 = rfs_1 - 1;
rfs = [num2str(rfs_1) ' : runs_left'];
[nothing, nothing] = unix(['echo "' rfs '" > multiple_runs']);
[nothing, nothing] = unix('rm -r CalcFold*');
[nothing, nothing] = unix('rm -r *.mat');
end
