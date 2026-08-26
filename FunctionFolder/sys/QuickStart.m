function Dummy = QuickStart(INPUT)
Names = fieldnames(INPUT);
Dummy = struct();
Num = size(Names, 1);
for i = 1:Num
Dummy = setfield(Dummy, Names{i}, '');
end
