function makeFigure_convexhull_301(generation, convex_hull, comp, points, points_new)
global ORG_STRUC
if size(ver('Octave'),1)
OctaveMode = 1;
else
OctaveMode = 0;
end
atomType = ORG_STRUC.atomType;
numIons  = ORG_STRUC.numIons;
x  = points(1,:);
y  = points(2,:);
x_ = points_new(1,:);
y_ = points_new(2,:);
if size(comp,2) == 2
A0=convex_hull(1,end-1); 
B0=convex_hull(2,end-1); 
end
if size(numIons,1) == 2
inConvexHull = zeros(2, size(convex_hull,1));
ConvexHull_atoms = zeros(size(convex_hull,1), size(atomType,2));
for i = 1:size(convex_hull,1)
N_atom  = convex_hull(i,1:end-2);
N_block = N_atom/numIons;
E = convex_hull(i,end-1);  
if size(numIons,2)==2
inConvexHull(1,i) = N_atom(2)/sum(N_atom); 
inConvexHull(2,i) = (E*sum(N_atom) - A0*sum(numIons(1,:))*N_block(1) ...
- B0*sum(numIons(2,:))*N_block(2)) / sum(N_atom); 
else
inConvexHull(1,i) = N_block(2)/sum(N_block); 
inConvexHull(2,i) = (E*sum(N_atom) - A0*sum(numIons(1,:))*N_block(1) ...
- B0*sum(numIons(2,:))*N_block(2)) / sum(N_block); 
end
ConvexHull_atoms(i,:) = N_atom;
end
try
h1 = figure;
set(gcf,'Visible','off');   
if OctaveMode == 0
scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y);
end
hold on;
box on;
[nothing, chRanking] = sort(inConvexHull(1,:)); 
line(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)), ...
'LineWidth', 1.0, 'Color', 'k');
if OctaveMode == 0
scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)),...
'MarkerEdgeColor','g','MarkerFaceColor','k');
else
scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)));
end
dy = (max(y)-min(y))/10;
xlim([min(x)-0.01, max(x)+0.01]);  
ylim([min(y)-dy,   max(y)+dy]);  
makeCaption(atomType, numIons);
set(gca,'FontSize',12,'linewidth',1.0);
print(h1,'-dpdf' ,        'extendedConvexHull.pdf');
makeLabel(atomType, ConvexHull_atoms, inConvexHull, dy/2);
if OctaveMode == 0
scatter(x_(1:end), y_(1:end), 'MarkerEdgeColor','k','MarkerFaceColor','r');
else
scatter(x_(1:end), y_(1:end));
end
print(h1,'-dpdf' ,       ['generation' num2str(generation)  '/extendedConvexHull.pdf']);
hold off
catch
end
end
function makeCaption(atomType, numIons)
if size(numIons,2)==2
name = ['Composition ratio: '];
name = [name megaDoof(atomType(2)) '/(' megaDoof(atomType(1)) '+' megaDoof(atomType(2)) ')'];
xlabel(name);
ylabel('Enthalpy of formation (eV/atom)');
else
name = '';
for type = 1:size(numIons,2)
num = numIons(1, type);
if num > 1
name = [name megaDoof(atomType(type)) num2str(num)];
elseif num>0
name = [name megaDoof(atomType(type))];
end
end
name = [name '+'];
for type = 1:size(numIons,2)
num = numIons(2, type);
if num > 1
name = [name megaDoof(atomType(type)) num2str(num)];
elseif num>0
name = [name megaDoof(atomType(type))];
end
end
xlabel(name);
ylabel('Enthalpy of formation (eV/block)');
end
function makeLabel(atomType, ConvexHull_atoms, ConvexHull, dy)
for loop1 = 1:size(ConvexHull_atoms, 1)
name = '';
for loop2 = 1:size(ConvexHull_atoms, 2)
if ConvexHull_atoms(loop1,loop2) > 0
label=['{' num2str(ConvexHull_atoms(loop1,loop2)) '}'];
name = [name megaDoof(atomType(loop2)) '_' label];
end
end
text(ConvexHull(1,loop1), ConvexHull(2,loop1)-dy, name, 'HorizontalAlignment', 'Center');
end
