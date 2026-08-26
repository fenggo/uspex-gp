function fitness = CalcFitness_310()
global POP_STRUC
global ORG_STRUC
global USPEX_STRUC
fitness = zeros(1,length(POP_STRUC.POPULATION));

disp('Calculating fitness ...');

for i = 1:length(POP_STRUC.POPULATION)
    if POP_STRUC.POPULATION(i).Enthalpies(end) < 9999
       if ORG_STRUC.optType == 1                              %%%%% entralpy
          factor = POP_STRUC.POPULATION(i).numMols/ORG_STRUC.numMols;
          fitness(i) = POP_STRUC.POPULATION(i).Enthalpies(end)/factor;
       elseif ORG_STRUC.optType == 2 
          fitness(i) = det(POP_STRUC.POPULATION(i).LATTICE)/sum(POP_STRUC.POPULATION(i).numIons);
       elseif ORG_STRUC.optType == 3                          %%%%%
          fitness(i) = -1*POP_STRUC.POPULATION(i).hardness;
       elseif ORG_STRUC.optType == 4                          %%%%%
          fitness(i) = -1*POP_STRUC.POPULATION(i).S_order;
       elseif ORG_STRUC.optType == 6 
          fitness(i) = -1*sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3;
       elseif ORG_STRUC.optType == 7 
          fitness(i) = -1*POP_STRUC.POPULATION(i).gap;
       elseif ORG_STRUC.optType == 8 
          Egc = 4; 
          if POP_STRUC.POPULATION(i).gap >= Egc
             fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^2; 
          else
             fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^6; 
          end
       elseif ORG_STRUC.optType == 9
          fitness(i) = -1*POP_STRUC.POPULATION(i).mag_moment;
       elseif ORG_STRUC.optType == 10 
          fitness(i) = -1*POP_STRUC.POPULATION(i).struc_entr;
       elseif ORG_STRUC.optType == 11                         %%%%%  density
          volume = det(POP_STRUC.POPULATION(i).LATTICE);
          density =  calcDensity( POP_STRUC.POPULATION(i).numIons, ORG_STRUC.atomType,volume);
          fitness(i) = -1*density;
       elseif ORG_STRUC.optType == 12                         %%%%%  density combined with entropy
          volume  = det(POP_STRUC.POPULATION(i).LATTICE);
          density =  calcDensity( POP_STRUC.POPULATION(i).numIons, ORG_STRUC.atomType,volume);
          factor  = (POP_STRUC.POPULATION(i).numMols/ORG_STRUC.numMols)*sum(ORG_STRUC.numMols);
          fitness(i) = POP_STRUC.POPULATION(i).Enthalpies(end)/factor - 0.5*density;
       elseif ORG_STRUC.optType == 13                         %%%%%  density combined with entropy
          volume  = det(POP_STRUC.POPULATION(i).LATTICE);
          density =  calcDensity( POP_STRUC.POPULATION(i).numIons, ORG_STRUC.atomType,volume);
          factor  = (POP_STRUC.POPULATION(i).numMols/ORG_STRUC.numMols)*sum(ORG_STRUC.numMols);
          fitness(i) = POP_STRUC.POPULATION(i).Enthalpies(end)/factor - density;
       elseif (ORG_STRUC.optType > 1100) & (ORG_STRUC.optType < 1112)  
          whichPara= mod(ORG_STRUC.optType,110);
          for i = 1 : length(POP_STRUC.POPULATION)
              if isempty(POP_STRUC.POPULATION(i).elasticProperties) | (POP_STRUC.POPULATION(i).elasticProperties(end)==0)
                 fitness(i)=NaN;
              else
                 fitness(i) = -1*POP_STRUC.POPULATION(i).elasticProperties(whichPara);
              end
          end
       end
    end
end

if ORG_STRUC.optType == 5 
   for i = 1 : length(POP_STRUC.POPULATION)-1
       for j = i+1 : length(POP_STRUC.POPULATION)
          dist_ij = cosineDistance(POP_STRUC.POPULATION(i).FINGERPRINT, POP_STRUC.POPULATION(j).FINGERPRINT, ORG_STRUC.weight);
          fitness(i) = fitness(i) + dist_ij^2;
          fitness(j) = fitness(j) + dist_ij^2;
       end
   end
   fitness = -sqrt(fitness);
end

fitness = ORG_STRUC.opt_sign*fitness; 
mean_fit = 0.0; 
n        = 0;
for i = 1 : length(fitness)                                %%%  fitness average
    if fitness(i) < 0
       mean_fit = mean_fit + fitness(i);
       n=n+1;
    end
end 
mean_fit = mean_fit/n;
variance = 0.0;                                           %%%  variance
for i = 1 : length(fitness)
    if fitness(i) < 0
       variance = variance + (fitness(i)-mean_fit)*(fitness(i)-mean_fit);
    end
end 
variance = variance/n;
tolerance = 3*variance;
if tolerance < 1.5
   tolerance = 1.5;
elseif tolerance > 2.0
    tolerance = 2.0;
end

for i = 1 : length(fitness)
    if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999
       fitness(i) = 100000;
    end
    if POP_STRUC.generation>1 && ORG_STRUC.optType == 12
       if fitness(i) < mean_fit - tolerance
          fitness(i) = 100001;
       end
    end 
end

% ===== GPML active learning with EI top-K acquisition =====
% Scan all valid crystals, compute EI, select top-K (K=3~5).
% If uncertainty > threshold, trigger DFT for all selected crystals
% in one batch → retrain GP → predict & update fitness for all valid.
% Parameters parsed from INPUT.txt: uspexkit gp --n= --data= --u=
if POP_STRUC.generation > 1
    valid_idx = find(fitness < 100000);
    if ~isempty(valid_idx)

        % --- parse parameters from INPUT.txt ---
        ncpu = 8;          % default
        dat_dir = 'data';   % default
        u_threshold = 0.04; % default
        top_k = 3;          % number of crystals to select per generation
        for ic = 1:length(ORG_STRUC.commandExecutable)
            cmd_str = ORG_STRUC.commandExecutable{ic};
            if ~isempty(strfind(cmd_str, 'uspexkit gp'))
                n_pos = strfind(cmd_str, '--n=');
                if ~isempty(n_pos)
                    rest = cmd_str(n_pos+4:end);
                    sp = find(rest == ' ', 1);
                    if isempty(sp), ncpu = str2num(rest);
                    else, ncpu = str2num(rest(1:sp-1)); end
                end
                d_pos = strfind(cmd_str, '--data=');
                if ~isempty(d_pos)
                    rest = cmd_str(d_pos+7:end);
                    sp = find(rest == ' ', 1);
                    if isempty(sp), dat_dir = rest;
                    else, dat_dir = rest(1:sp-1); end
                end
                u_pos = strfind(cmd_str, '--u=');
                if ~isempty(u_pos)
                    rest = cmd_str(u_pos+4:end);
                    sp = find(rest == ' ', 1);
                    if isempty(sp), u_threshold = str2num(rest);
                    else, u_threshold = str2num(rest(1:sp-1)); end
                end
                break;
            end
        end

        % --- EI top-K selection: scan all valid crystals ---
        gp_file = [ORG_STRUC.homePath '/CalcFold1/gp.csv'];
        if exist(gp_file, 'file')
            data = csvread(gp_file, 1, 0);    % skip header
            % gp.csv: csvread 跳过 header 后
            % data 列: 1=bodyCount(自身ID), 2=最近邻索引, 3=residual, 4=den_min,
            %   5=den_rf, 6=den_gp, 7=uncertainty, 8=eng_min, 9=eng_pred, 10=uncert_eng
            idx_col = data(:, 1);     % bodyCount = 晶体自身ID
            den_gp_col = data(:, 6);  % density_gp
            uncert_col = data(:, 7);  % uncertainty
            resid_col = data(:, 3);   % residual

            % EI = Expected Improvement (无超参数, 自动平衡探索)
            % gp.csv col1 = bodyCount, col2 = 最近邻索引, col3 = residual

            % Step 0: compute f_best = max density_gp among valid crystals
            f_best = -inf;
            filtered_count = 0;
            for vi = 1:length(valid_idx)
                cnum = POP_STRUC.POPULATION(valid_idx(vi)).Number;
                match = find(idx_col == cnum);
                if ~isempty(match)
                    r = match(1);
                    if resid_col(r) > 10
                        filtered_count = filtered_count + 1;
                        continue;
                    end
                    if den_gp_col(r) > f_best
                        f_best = den_gp_col(r);
                    end
                end
            end
            if filtered_count > 0
                fprintf('  Filtered out %d crystals with residual > 10\n', filtered_count);
            end

            % Step 1: compute EI for all valid crystals, collect into arrays
            ei_list = [];
            ei_crystal_list = [];
            ei_uncert_list = [];
            ei_den_gp_list = [];
            for vi = 1:length(valid_idx)
                cnum = POP_STRUC.POPULATION(valid_idx(vi)).Number;
                match = find(idx_col == cnum);
                if ~isempty(match)
                    r = match(1);
                    if resid_col(r) > 10
                        continue;
                    end
                    mu = den_gp_col(r);
                    sigma = uncert_col(r);

                    % EI formula for maximization:
                    %   Z = (mu - f_best) / sigma
                    %   EI = (mu - f_best) * Phi(Z) + sigma * phi(Z)
                    diff = mu - f_best;
                    if sigma < 1e-10
                        ei = max(0, diff);
                    else
                        z = diff / sigma;
                        Phi_z = 0.5 * erfc(-z / sqrt(2));
                        phi_z = exp(-0.5 * z^2) / sqrt(2 * pi);
                        ei = diff * Phi_z + sigma * phi_z;
                    end

                    ei_list = [ei_list; ei];
                    ei_crystal_list = [ei_crystal_list; cnum];
                    ei_uncert_list = [ei_uncert_list; sigma];
                    ei_den_gp_list = [ei_den_gp_list; mu];
                end
            end

            % Step 2: sort by EI descending, select top-K with uncertainty > threshold
            [~, sort_idx] = sort(ei_list, 'descend');
            n_sel = min(top_k, length(ei_list));
            selected_crystals = ei_crystal_list(sort_idx(1:n_sel));
            selected_uncerts = ei_uncert_list(sort_idx(1:n_sel));
            selected_den_gps = ei_den_gp_list(sort_idx(1:n_sel));
            selected_eis = ei_list(sort_idx(1:n_sel));

            % Only keep crystals with uncertainty > threshold
            select_mask = selected_uncerts > u_threshold;
            selected_crystals = selected_crystals(select_mask);
            selected_uncerts = selected_uncerts(select_mask);
            selected_den_gps = selected_den_gps(select_mask);
            selected_eis = selected_eis(select_mask);

            n_dft = length(selected_crystals);
            if n_dft > 0
                fprintf('EI top-%d selection: %d/%d crystals above threshold (%.4f)\n', ...
                    top_k, n_dft, n_sel, u_threshold);
                for s = 1:n_dft
                    fprintf('  %d: crystal %d  EI=%.6f  den_gp=%.4f  sigma=%.4f  f_best=%.4f\n', ...
                        s, selected_crystals(s), selected_eis(s), ...
                        selected_den_gps(s), selected_uncerts(s), f_best);
                end

                % Step 3: generate trajectory (once for all)
                cmd = ['cd ' ORG_STRUC.resFolder ' && uspexkit traj'];
                fprintf('  Running: %s\n', cmd);
                system(cmd);

                % Step 4: DFT calculation for all selected crystals
                ids_str = '';
                for s = 1:n_dft
                    if s == 1
                        ids_str = num2str(selected_crystals(s));
                    else
                        ids_str = [ids_str ' ' num2str(selected_crystals(s))];
                    end
                end
                cmd = sprintf('cd %s && uspexkit calc --ncpu=%d --ids=''%s'' --dat=%s', ...
                              ORG_STRUC.resFolder, ncpu, ids_str, dat_dir);
                fprintf('  Running: %s\n', cmd);
                system(cmd);

                % Step 5: GP model retrained → pred valid crystals (residual <= 10)
                pred_ids_str = '';
                pred_id_count = 0;
                for vi = 1:length(valid_idx)
                    cnum = POP_STRUC.POPULATION(valid_idx(vi)).Number;
                    match = find(idx_col == cnum);
                    if ~isempty(match) && resid_col(match(1)) <= 10
                        if pred_id_count == 0
                            pred_ids_str = num2str(cnum);
                        else
                            pred_ids_str = [pred_ids_str ' ' num2str(cnum)];
                        end
                        pred_id_count = pred_id_count + 1;
                    end
                end

                if pred_id_count > 0
                    pred_log = [ORG_STRUC.resFolder '/density_predict.log'];
                    if exist(pred_log, 'file')
                        delete(pred_log);
                    end
                    cmd = sprintf('cd %s && uspexkit pred --step=100 --ncpu=%d --den=0.0 --ids=''%s'' --dat=%s', ...
                                  ORG_STRUC.resFolder, ncpu, pred_ids_str, dat_dir);
                    fprintf('  Running: %s\n', cmd);
                    system(cmd);

                    % read pred results, update fitness (skip residual > 10)
                    if exist(pred_log, 'file')
                        % density_predict.log 是空格分隔，非 CSV
                        pred_raw = importdata(pred_log);
                        pred_data = pred_raw.data;  % numeric matrix
                        pred_id_col = pred_data(:, 1);
                        pred_density_col = pred_data(:, 5);

                        updated_count = 0;
                        skipped_count = 0;
                        for vi = 1:length(valid_idx)
                            idx = valid_idx(vi);
                            cnum = POP_STRUC.POPULATION(idx).Number;
                            match = find(idx_col == cnum);
                            if ~isempty(match) && resid_col(match(1)) > 10
                                skipped_count = skipped_count + 1;
                                continue;
                            end
                            pred_match = find(pred_id_col == cnum);
                            if ~isempty(pred_match)
                                fitness(idx) = -1 * pred_density_col(pred_match(1));
                                updated_count = updated_count + 1;
                            end
                        end
                        fprintf('  Updated fitness for %d/%d valid crystals (skipped %d with residual > 10)\n', ...
                            updated_count, length(valid_idx), skipped_count);
                    end
                else
                    fprintf('  No valid crystals with residual <= 10, skipping pred\n');
                end
            else
                fprintf('  -> No crystals above uncertainty threshold (%.4f), skipping DFT and pred\n', u_threshold);
            end
        end
    end
end