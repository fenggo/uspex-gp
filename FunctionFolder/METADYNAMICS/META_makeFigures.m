function META_makeFigures(Nsteps)
global USPEX_STRUC
if size(ver('Octave'),1)
OctaveMode = 1;
else
OctaveMode = 0;
end
N             = length(USPEX_STRUC.POPULATION);
N_gen         = length(USPEX_STRUC.GENERATION);
enth          = zeros(N,Nsteps);
BestE         = zeros(N_gen,1);
BestE_relaxed = zeros(N_gen,1);
lat_basic     = zeros(N_gen,6);
for i=1:N
numIons   = USPEX_STRUC.POPULATION(i).numIons;
enth(i,:) = USPEX_STRUC.POPULATION(i).Enthalpies(1:Nsteps)/sum(numIons);
end
for i=1:N_gen
ID(i) = i-1;
BestE(i)         = USPEX_STRUC.GENERATION(i).Best_enth;
BestE_relaxed(i) = USPEX_STRUC.GENERATION(i).Best_enth_relaxed;
lat_basic(i,:)   = USPEX_STRUC.GENERATION(i).lat_basic;
end
try
[nothing, nothing] = unix(['rm BestEnthalpy.pdf']);
x = ID;
y1 = BestE;
y2 = BestE_relaxed;
h = figure;
set(gcf,'Visible','off');   
plot(x,y2,x,y1, '--rs','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',5);
xlabel('Generation number');
ylabel('Best enthalpy (eV/atom)');
print(h,'-dpdf', 'BestEnthalpy.pdf');
catch
end
try
h = figure;
set(gcf,'Visible','off');   
e_complete = enth;
sub = ceil((Nsteps-1)/2);
[nothing, nothing] = unix(['rm E_series.pdf']);
for g = 1 : Nsteps - 1
x = []; y = [];
x = e_complete(:,g);
y = e_complete(:,g+1);
subplot(sub,2,g);
if OctaveMode == 0
scatter(x,y,10,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y,10);
end
xlabel(['E' num2str(g)]);
ylabel(['E' num2str(g+1)]);
end
print(h,'-dpdf', 'E_series.pdf');
catch
end
y = lat_basic;
x = ID;
h = figure;
set(gcf,'Visible','off');   
y_tmp = y(:,1);
subplot(3,2,1);
scatter_size = 15;
if OctaveMode == 0
scatter(x,y_tmp,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y_tmp, scatter_size);
end
xlabel('N');
ylabel('a_{11}');
y_tmp = y(:,2);
subplot(3,2,3);
if OctaveMode == 0
scatter(x,y_tmp,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y_tmp, scatter_size);
end
xlabel('N');
ylabel('a_{22}');
y_tmp = y(:,3);
subplot(3,2,5);
if OctaveMode == 0
scatter(x,y_tmp,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y_tmp, scatter_size);
end
xlabel('N');
ylabel('a_{33}');
y_tmp = y(:,4);
subplot(3,2,2);
if OctaveMode == 0
scatter(x,y_tmp,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y_tmp, scatter_size);
end
xlabel('N');
ylabel('a_{21}');
y_tmp = y(:,5);
subplot(3,2,4);
if OctaveMode == 0
scatter(x,y_tmp,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y_tmp, scatter_size);
end
xlabel('N');
ylabel('a_{31}');
y_tmp = y(:,6);
subplot(3,2,6);
if OctaveMode == 0
scatter(x,y_tmp,'MarkerEdgeColor','k','MarkerFaceColor','g');
else
scatter(x,y_tmp, scatter_size);
end
xlabel('N');
ylabel('a_{32}');
print(h,'-dpdf', 'Lattice_vs_N.pdf');
