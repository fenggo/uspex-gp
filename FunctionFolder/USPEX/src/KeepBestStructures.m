function addon_diff = KeepBestStructures()
global ORG_STRUC
global POP_STRUC
global OFF_STRUC

keepBestHM = ORG_STRUC.keepBestHM;
decentRank = min(round(length(POP_STRUC.ranking)*ORG_STRUC.bestFrac), length(POP_STRUC.ranking)-POP_STRUC.bad_rank );
if keepBestHM > decentRank
    keepBestHM = decentRank;
end

fitness = POP_STRUC.fitness;
first_mean = mean(fitness(POP_STRUC.ranking(1:round(ORG_STRUC.bestFrac*length(POP_STRUC.ranking)))));

% SOAP 混合距离参数
soapAlpha = 0.9;
if isfield(ORG_STRUC, 'soapAlpha')
    soapAlpha = ORG_STRUC.soapAlpha;
end

if ~isempty(POP_STRUC.convex_hull)
    for i = 1 : size(POP_STRUC.convex_hull,1)
        IND = POP_STRUC.convex_hull(i,end);
        if IND > 0
            if isempty(OFF_STRUC.POPULATION(end).numIons)
                OFF_STRUC.POPULATION(end) = POP_STRUC.POPULATION(IND);
            else
                OFF_STRUC.POPULATION(end+1) = POP_STRUC.POPULATION(IND);
            end
            OFF_STRUC.POPULATION(end).Parents = [];
            if ORG_STRUC.reoptOld
                OFF_STRUC.POPULATION(end).Step = max(length(ORG_STRUC.abinitioCode)-1,1);
            else
                OFF_STRUC.POPULATION(end).Step = length(ORG_STRUC.abinitioCode)+1;
            end
            info_parents = struct('parent',{}, 'enthalpy',{});
            info_parents(1).parent = num2str(POP_STRUC.POPULATION(IND).Number);
            info_parents.enthalpy = POP_STRUC.POPULATION(IND).Enthalpies(end)/sum(POP_STRUC.POPULATION(IND).numIons);
            OFF_STRUC.POPULATION(end).Parents = info_parents;
            OFF_STRUC.POPULATION(end).howCome = 'convexHull';
        end
    end
    if ORG_STRUC.dimension==3
        keepBestHM = 0;
    else
        keepBestHM = keepBestHM - size(POP_STRUC.convex_hull,1);
    end
    if keepBestHM < 0
        keepBestHM = 0;
    end
end

if (ORG_STRUC.softMutOnly(POP_STRUC.generation) == 1)
    for i = 1 : length(POP_STRUC.ranking) - POP_STRUC.bad_rank
        IND = POP_STRUC.ranking(i);
        OFF_STRUC.POPULATION(end+1) = POP_STRUC.POPULATION(IND);
        OFF_STRUC.POPULATION(end).Parents = [];
        info_parents = struct('parent',{}, 'enthalpy', {});
        info_parents(1).parent = num2str(POP_STRUC.POPULATION(IND).Number);
        info_parents.enthalpy = POP_STRUC.POPULATION(IND).Enthalpies(end)/sum(POP_STRUC.POPULATION(IND).numIons);
        OFF_STRUC.POPULATION(end).Parents = info_parents;
        OFF_STRUC.POPULATION(end).howCome = 'keptBest';
        if ORG_STRUC.reoptOld
            OFF_STRUC.POPULATION(end).Step = max(length(ORG_STRUC.abinitioCode)-1,1);
        else
            OFF_STRUC.POPULATION(end).Step = length(ORG_STRUC.abinitioCode)+1;
        end
    end
    keepBestHM = 0;
    addon_diff = length(POP_STRUC.ranking) - POP_STRUC.bad_rank;
    disp([num2str(i) ' structures were kept from the convex hull']);
else
    tolerance = ORG_STRUC.toleranceBestHM;
    deltaTol = tolerance/2;
    doneOr = 0;
    if ORG_STRUC.dynamicalBestHM == 2
        tolerance = 0.5;
        deltaTol = tolerance/2;
    end
    if (ORG_STRUC.softMutOnly(POP_STRUC.generation) == 1)
        doneOr = 1;
        chosen = zeros(1,length(POP_STRUC.ranking));
    end
    while doneOr == 0
        chosen = zeros(1,length(POP_STRUC.ranking));
        addon = 1;
        addon_diff = 0;
        while addon_diff < keepBestHM
            if addon > length(POP_STRUC.ranking)
                break;
            end
            if (addon > decentRank) & (ORG_STRUC.dynamicalBestHM > 0)
                break;
            end
            if ~isempty(POP_STRUC.convex_hull)
                if ~isempty(find(POP_STRUC.convex_hull(:,end) == POP_STRUC.ranking(addon)))
                    addon = addon + 1;
                    continue;
                end
            end
            good_structure = 1;
            if (ORG_STRUC.doFing) & (POP_STRUC.generation > ORG_STRUC.doFing-1)
                % --- 检查与已保留结构的相似度 ---
                if ~isempty(POP_STRUC.convex_hull)
                    for j = 1 : size(POP_STRUC.convex_hull,1)
                        fing = OFF_STRUC.POPULATION(end-j+1).FINGERPRINT;
                        if ~isempty(fing)
                            soapA = []; soapB = [];
                            if isfield(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)), 'soap_fp')
                                soapA = POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).soap_fp;
                            end
                            if isfield(OFF_STRUC.POPULATION(end-j+1), 'soap_fp')
                                soapB = OFF_STRUC.POPULATION(end-j+1).soap_fp;
                            end
                            if isempty(soapA) || isempty(soapB)
                                cos_dist = cosineDistance(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).FINGERPRINT, fing, ORG_STRUC.weight);
                            else
                                cos_dist = hybridCosineDistance(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).FINGERPRINT, fing, soapA, soapB, ORG_STRUC.weight, soapAlpha);
                            end
                            if (cos_dist < tolerance)
                                if (ORG_STRUC.varcomp == 1) && (ORG_STRUC.dimension == 0)
                                    numions1 = OFF_STRUC.POPULATION(end-j+1).numIons;
                                    numions2 = POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).numIons;
                                    if numions1==numions2
                                        good_structure = 0;
                                    end
                                else
                                    good_structure = 0;
                                end
                            end
                        end
                    end
                end
                % --- 检查与本轮已选结构的相似度 ---
                for j = 1 : length(POP_STRUC.ranking)
                    if chosen(j) == 1
                        soapA2 = []; soapB2 = [];
                        if isfield(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)), 'soap_fp')
                            soapA2 = POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).soap_fp;
                        end
                        if isfield(POP_STRUC.POPULATION(POP_STRUC.ranking(j)), 'soap_fp')
                            soapB2 = POP_STRUC.POPULATION(POP_STRUC.ranking(j)).soap_fp;
                        end
                        if isempty(soapA2) || isempty(soapB2)
                            cos_dist = cosineDistance(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).FINGERPRINT, POP_STRUC.POPULATION(POP_STRUC.ranking(j)).FINGERPRINT, ORG_STRUC.weight);
                        else
                            cos_dist = hybridCosineDistance(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).FINGERPRINT, POP_STRUC.POPULATION(POP_STRUC.ranking(j)).FINGERPRINT, soapA2, soapB2, ORG_STRUC.weight, soapAlpha);
                        end
                        if (cos_dist < tolerance)
                            if (ORG_STRUC.varcomp == 1) && (ORG_STRUC.dimension == 0)
                                numions1 = POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).numIons;
                                numions2 = POP_STRUC.POPULATION(POP_STRUC.ranking(j)).numIons;
                                if numions1==numions2
                                    good_structure = 0;
                                end
                            else
                                good_structure = 0;
                            end
                        end
                    end
                end
            end
            if ((length(POP_STRUC.POPULATION)-addon)-(keepBestHM-addon_diff) <= 1) & (ORG_STRUC.dynamicalBestHM == 0)
                good_structure = 1;
            end
            if ORG_STRUC.dynamicalBestHM
                if addon > 1
                    if fitness(POP_STRUC.ranking(addon)) > first_mean + power(var(fitness(POP_STRUC.ranking(1:round(ORG_STRUC.bestFrac*length(POP_STRUC.ranking))))), 0.5)
                        break;
                    end;
                end;
            end
            if good_structure
                chosen(addon) = 1;
                addon_diff = addon_diff + 1;
            end
            addon = addon+1;
        end
        if ORG_STRUC.dynamicalBestHM < 2
            doneOr = 1;
        else
            if sum(chosen(:)) == keepBestHM
                doneOr = 1;
            elseif sum(chosen(:)) < keepBestHM
                tolerance = tolerance - deltaTol;
            else
                tolerance = tolerance + deltaTol;
            end
            deltaTol = deltaTol/2;
            if deltaTol < 0.000001
                doneOr = 1;
            end
        end
    end
    for addon = 1 : length(POP_STRUC.ranking)
        if chosen(addon)
            if isempty(OFF_STRUC.POPULATION(end).numIons)
                OFF_STRUC.POPULATION(end)   = POP_STRUC.POPULATION(POP_STRUC.ranking(addon));
            else
                OFF_STRUC.POPULATION(end+1) = POP_STRUC.POPULATION(POP_STRUC.ranking(addon));
            end
            OFF_STRUC.POPULATION(end).Parents = [];
            info_parents = struct('parent',{}, 'enthalpy', {});
            info_parents(1).parent = num2str(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).Number);
            info_parents.enthalpy = POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).Enthalpies(end);
            if strcmp(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).numIons, 'N/A') ~= 1
                info_parents.enthalpy = info_parents.enthalpy/sum(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).numIons);
            end
            OFF_STRUC.POPULATION(end).Parents = info_parents;
            OFF_STRUC.POPULATION(end).howCome = 'keptBest';
            if ORG_STRUC.reoptOld
                OFF_STRUC.POPULATION(end).Step = max(length(ORG_STRUC.abinitioCode)-1,1);
            else
                OFF_STRUC.POPULATION(end).Step = length(ORG_STRUC.abinitioCode)+1;
            end
        end
    end
    disp([num2str(addon_diff) ' structures were kept as best']);
end
