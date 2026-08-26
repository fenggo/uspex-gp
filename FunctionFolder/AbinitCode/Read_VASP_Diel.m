function [d_s] = Read_VASP_Diel()
d_s = zeros(1,6);
d_s_ion = zeros(1,6);
[fileHandle, msg] = fopen('OUTCAR');
while 1
tline = fgetl(fileHandle);
if ~ischar(tline)
break
end
if ~isempty(strfind(tline, 'MACROSCOPIC STATIC DIELECTRIC TENSOR ('))
tline = fgetl(fileHandle);
tline1 = fgetl(fileHandle);
tline2 = fgetl(fileHandle);
tline3 = fgetl(fileHandle);
e_tensor = str2num(vertcat(tline1, tline2, tline3));
d_s(1) = e_tensor(1,1);
d_s(2) = e_tensor(2,2);
d_s(3) = e_tensor(3,3);
d_s(4) = e_tensor(1,2);
d_s(5) = e_tensor(1,3);
d_s(6) = e_tensor(2,3);
end
if ~isempty(strfind(tline, 'MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC'))
tline = fgetl(fileHandle);
tline1 = fgetl(fileHandle);
tline2 = fgetl(fileHandle);
tline3 = fgetl(fileHandle);
e_tensor = str2num(vertcat(tline1, tline2, tline3));
d_s_ion(1) = e_tensor(1,1);
d_s_ion(2) = e_tensor(2,2);
d_s_ion(3) = e_tensor(3,3);
d_s_ion(4) = e_tensor(1,2);
d_s_ion(5) = e_tensor(1,3);
d_s_ion(6) = e_tensor(2,3);
end
end
status = fclose(fileHandle);
d_s = d_s + d_s_ion;
[tmp, imaginary_f] = unix('./getStuff OUTCAR f/i 4'); 
if ~isempty(imaginary_f)
i_f = abs(str2num(imaginary_f));
if (i_f(end) > 0.3) | (i_f(end-1) > 0.3) | (i_f(end-2) > 0.3) 
d_s = zeros(1,6);
end
end
