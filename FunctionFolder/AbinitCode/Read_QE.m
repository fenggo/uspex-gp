function target = Read_QE(flag, ID)
Ry = 13.6056923; 
bohr= 0.529177249;
if     flag == -1  
[nothing, results] = unix('grep "JOB DONE" output');
if isempty(results)
disp('Quantum Espresso is not completely Done');
[nothing, nothing] = unix(['cp output ERROR-OUTPUT-' ID]);
target = 0;
else
target = 1;
end
elseif flag == 0
[nothing, results1] = unix('./getStuff output "enthalpy new" 5'); 
[nothing, results2] = unix('./getStuff output "!    total energy" 5'); 
if isempty(results1) & isempty(results2)
disp('Quantum Espresso 1st SCF is not done');
[nothing, nothing] = unix(['cp output ERROR-OUTPUT-' ID]);
target = 0;
else
target = 1;
end
elseif flag ==  1 
[nothing, results1] = unix('./getStuff output "enthalpy new" 5 |tail -1'); 
if ~isempty(results1) 
target = str2num(results1);
else
[nothing, results2] = unix('./getStuff output "!    total energy" 5 |tail -1'); 
target = str2num(results2);
end
target = target*Ry;
elseif flag ==  2 
target = callAWK('QE_pres.awk','output');
elseif flag ==  7
try
[nothing, numStr] = unix(['grep "atoms/cell" output |cut -d"=" -f2']);
numIons= sum( str2num(numStr) );
force_orig = callAWK('QE_force.awk', 'output', ['num=', num2str(numIons)]);
if isempty(force_orig)
target = [];
else
target = force_orig*Ry/bohr;
end
catch
target = [];
end
end
