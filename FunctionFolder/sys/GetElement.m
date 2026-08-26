function atomType1 = GetElement(N_type, atomType)
atomType1 = zeros(1,N_type);
c1 = findstr(atomType, ' ');
c = sort(str2num(['0 ' num2str(c1)]));
c(end+1) = length(atomType) + 1;
ind1 = 1;
for i = 2 : length(c)
if c(i-1)+1 > c(i)-1
continue
end
tmp = atomType(c(i-1)+1 : c(i)-1);
if ~isempty(str2num(tmp)) 
atomType1(ind1) = str2num(tmp);
else
atomType1(ind1)=elementID(tmp);  
end
ind1 = ind1 + 1;
end
