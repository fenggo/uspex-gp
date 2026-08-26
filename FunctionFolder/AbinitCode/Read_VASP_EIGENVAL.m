function [eGap] = read_EIGENVAL
try
handle = fopen('EIGENVAL');
for i = 1 : 6
tmp = fgetl(handle);
end
l6 = str2num(tmp);
eN = l6(1);
filledBand = floor(eN/2); 
nK = l6(2); 
nB = l6(3); 
for i = 1 : nK
tmp = fgetl(handle); 
tmp = fgetl(handle); 
bands = fscanf(handle, '%g %g', [2 nB]);
bands = bands';
if i == 1
gapStart = bands(filledBand, 2);
gapEnd = bands(filledBand+1, 2);
else
if gapStart < bands(filledBand, 2);
gapStart = bands(filledBand, 2);
end
if gapEnd > bands(filledBand+1, 2);
gapEnd = bands(filledBand+1, 2);
end
end
tmp = fgetl(handle); 
end
fclose(handle);
eGap = gapEnd - gapStart;
if eGap < 0
eGap = 0; 
end
catch
eGap = 0;
end