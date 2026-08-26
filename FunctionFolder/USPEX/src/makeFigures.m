function makeFigures(pickUpNCount, Nsteps, ConstLat)
global USPEX_STRUC

if size(ver('Octave'),1)
   OctaveMode = 1;
else
   OctaveMode = 0;
end

N         = length(USPEX_STRUC.POPULATION);
enth      = zeros(N,Nsteps);
vol       = zeros(N,1);
fit       = zeros(N,1);
ID        = zeros(N,1);

for i=1:N
    numIons   = USPEX_STRUC.POPULATION(i).numIons;
    fit(i)    = USPEX_STRUC.POPULATION(i).Fitness;
    enth(i,:) = USPEX_STRUC.POPULATION(i).Enthalpies/sum(numIons);
    vol(i) =det(USPEX_STRUC.POPULATION(i).LATTICE)/sum(numIons);
    ID(i) = i;
end

try
   [a,b]=unix(['rm Energy_vs_N.pdf']);
   x = []; y = [];
   x = ID;
   y = enth(:,Nsteps);
   h = figure;
   set(gcf,'Visible','off');   
   if OctaveMode == 0
      scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
   else
      scatter(x,y);
   end
   xlabel('Structure number');
   ylabel('Enthalpy');
   print(h,'-dpdf', 'Energy_vs_N.pdf');
   catch
end

try
   [a,b]=unix(['rm Fitness_vs_N.pdf']);
   x = []; y = [];
   x = ID;
   y = fit;
   y_max = max(y(find(y<1000))); 
   y(find(y>=10000))=y_max;
   h = figure;
   set(gcf,'Visible','off');   
   if OctaveMode == 0
      scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
   else
      scatter(x,y);
   end
   xlabel('Structure number');
   ylabel('Fitness');
   print(h,'-dpdf','-r120','Fitness_vs_N.pdf');
   catch
end

if ConstLat~=1
   try
      [a,b]=unix(['rm Energy_vs_Volume.pdf']);
      x = []; y = [];
      x = enth(:,Nsteps);
      y = vol;
      h = figure;
      set(gcf,'Visible','off');   
      if OctaveMode == 0
         scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
      else
         scatter(x,y);
      end
      xlabel('Enthalpy');
      ylabel('Volume');
      print(h,'-dpdf','Energy_vs_Volume.pdf');
      catch
   end
end

try
   e_complete = enth;
   sub = ceil((Nsteps-1)/2);
   [a,b]=unix(['rm E_series.pdf']);
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

try
   [a, b]=unix(['rm Variation-Operators.pdf']);
   h = figure;
   set(gcf,'Visible','off');   
   N_count = N-pickUpNCount;
   x = ID(N-N_count+1:N);
   y = zeros(N_count,1);
   s = zeros(N_count,1);  
   c = zeros(N_count,1);  
   handle = fopen('origin');
   tmp = fgetl(handle);
   for i=1:N_count
       tmp = fgetl(handle);
       if tmp ~= -1 && isempty(findstr(tmp, 'Random'))
          if ~isempty(findstr(tmp, 'LatMutated'))
             c(i) = 0;
             s(i) = 30;
          elseif ~isempty(findstr(tmp, 'CoorMutate'))  
             c(i) = 1;
             s(i) = 30;
          elseif ~isempty(findstr(tmp, 'softmutate'))
             c(i) = 1;
             s(i) = 30;
          elseif ~isempty(findstr(tmp, 'Rotated'))
             c(i) = 2;
             s(i) = 30;
          elseif ~isempty(findstr(tmp, 'Permutate'))
             c(i) = 3;
             s(i) = 30;
          elseif ~isempty(findstr(tmp, 'Heredity'))
             c(i) = 4;
             s(i) = 30;
          end
             CE = str2num(tmp(17:26));
             PE = str2num(tmp(27:36));
             y(i) = CE-PE;
       end
   end
   status = fclose(handle);
   subplot(2,1,1);
   NonZeroLocations = s ~= 0;
   
   if OctaveMode == 0
      scatter(x(NonZeroLocations), y(NonZeroLocations), s(NonZeroLocations), c(NonZeroLocations),'filled','MarkerEdgeColor','k');
   else
      scatter(x(NonZeroLocations), y(NonZeroLocations), s(NonZeroLocations), c(NonZeroLocations));
   end
   
   xlabel('Structure index');
   ylabel('E_{Child} - E_{parent}');
   hcb = colorbar('YTick',[0:1:4],'YTickLabel',{'Latmut.','AtomMut.','Rotation','Permut.','Heredity'});
   set(hcb,'YTickMode','manual')
   [a, b] = unix(['grep frac OUTPUT.txt |cut -c60-65']);
   
   if ~isempty(b)
      str = {'Heredity', 'Random', 'Softmutation', 'Permutation', 'Latmutation', 'Rotation', 'Transmutation'};
      num_of_operators = size(str, 2);
      x=[]; y=[];
      frac = str2num(b);
      N = size(frac,1)/num_of_operators;
      x = 1:1:N;
      y = zeros(num_of_operators, N);
      for i = 1:N
          for j = 1:num_of_operators
              y(j,i) = frac(num_of_operators * (i - 1) + j);
          end
      end
      subplot(2,1,2);
      plot(x,y(1,:),x,y(2,:),x,y(3,:),x,y(4,:),x,y(5,:),x,y(6,:),x,y(7,:));
      legend(str);
      xlabel('Generation number');
      ylabel('Fraction of each variation');
      if N>1
         axis([1 N 0 1]);
      end
   end
   print(h,'-dpdf', 'Variation-Operators.pdf');
   catch
end
