function target = Read_CASTEP(flag)
if flag == 0
[nothing, results] = unix('grep ''Final Enthalpy'' cstp.castep | awk -F= ''{print $2}'' | awk ''{print $1}'' ');
if isempty(results) | isempty(str2num(results))
disp('CASTEP relaxation is not done');
target = 0;
else
target = 1;
end
elseif flag ==  1 
[nothing, results] = unix('grep ''Final Enthalpy'' cstp.castep | awk -F= ''{print $2}'' | awk ''{print $1}'' ');
fitnesses = str2num(results);
target    = fitnesses(end);
end
