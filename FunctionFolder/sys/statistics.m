function statistics(target_enthalpy)
function dirs = getAllDirs(path, dirs)
listing = dir(path);
for i = 1:size(listing,1)
if listing(i).isdir & strcmp(listing(i).name, '..') == 0 & ...
strcmp(listing(i).name, '.') == 0
dirs{end+1} = [path '/' listing(i).name];
dirs = getAllDirs([path '/' listing(i).name], dirs);
end
end
end
dirs = getAllDirs(pwd, {});
results_dirs = {};
for i = 1:size(dirs,2)
dirs{i};
if strfind(dirs{i}, 'results') > 0
[pathstr, name, ext] = fileparts(dirs{i});
if strfind(name, 'results') > 0
%disp(dirs{i});
results_dirs{end+1} = dirs{i};
end
end
end
USPEX_mat = {};
for i = 1:size(results_dirs,2)
if exist([results_dirs{i} '/' 'USPEX.mat'], 'file')
%disp([results_dirs{i} '/' 'USPEX.mat']);
USPEX_mat{end+1} = [results_dirs{i} '/' 'USPEX.mat'];
end
end
number_of_mat = size(USPEX_mat,2);
disp(' ');
disp(['Number of files to be processed: ' num2str(number_of_mat)]);
disp(['Target enthalpy: ' num2str(target_enthalpy)]);
disp(' ');
if number_of_mat > 0
found_generations = -1*ones(number_of_mat,1);
found_populations = -1*ones(number_of_mat,1);
found_enthalpies  = nan(number_of_mat,1);
for i = 1:number_of_mat
%disp(USPEX_mat{i})
global USPEX_STRUC
load(USPEX_mat{i});
pop_enth_full = nan(size(USPEX_STRUC.POPULATION(:)));
pop_enth      = [];
for j=1:size(USPEX_STRUC.POPULATION(:))
howcome = USPEX_STRUC.POPULATION(j).howCome;
if ~strcmp(howcome, 'keptBest') & ~strcmp(howcome, 'convexHull')
pop_enth_full(j) = USPEX_STRUC.POPULATION(j).Enthalpies(end);
pop_enth         = [pop_enth; USPEX_STRUC.POPULATION(j).Enthalpies(end)];
end
end
found_id_full = find(pop_enth_full <= target_enthalpy);
found_id      = find(pop_enth      <= target_enthalpy);
if ~isempty(found_id)
first_id_full = found_id_full(1);
first_id      = found_id(1);
first_gen  = USPEX_STRUC.POPULATION(first_id_full).gen;
first_enth = pop_enth(first_id);
found_generations(i) = first_gen;
found_populations(i) = first_id;
found_enthalpies(i)  = first_enth;
str = sprintf('Generation: %2i    Number: %4i    Enthalpy: %10.4f    Mat-file: %s', ...
first_gen,         first_id,      first_enth,         USPEX_mat{i});
disp(str);
end
clear('USPEX_STRUC');
end
disp(' ');
% disp([found_populations    found_generations    found_enthalpies]);
data_pop = [];
data_gen = [];
format   = [];
for i=1:number_of_mat
format   = [format '
data_pop = [data_pop found_populations(i)];
data_gen = [data_gen found_generations(i)];
end
str_pop = sprintf(['Found structures numbers : ' format], data_pop);
str_gen = sprintf(['Found generations numbers: ' format], data_gen);
disp(str_pop);
disp(str_gen)
positive_gen = find(found_generations > 0);
positive_pop = find(found_populations > 0);
success_rate = size(positive_gen,1)/number_of_mat*100;
mean_gens = ceil(mean(found_generations(positive_gen)));
mean_pops = ceil(mean(found_populations(positive_pop)));
std_pops  = ceil( std(found_populations(positive_pop)));
disp(' ');
disp(['Success rate: ' num2str(success_rate) ' percent']);
disp(['Average number of generations to get E=' num2str(target_enthalpy) ': ' num2str(mean_gens)]);
disp(['Average number of structures  to get E=' num2str(target_enthalpy) ': ' num2str(mean_pops)]);
disp(['Standard deviation: ' num2str(std_pops)]);
disp(' ');
end
end
