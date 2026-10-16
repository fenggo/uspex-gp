function ReadJobs_310()
global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    if ~isempty (whichInd)
        Step = POP_STRUC.POPULATION(whichInd).Step;
        disp(['Structure' num2str(whichInd) ' step' num2str(Step) ' at CalcFold' num2str(indic) ]);
        if POP_STRUC.POPULATION(whichInd).JobID
            if (ORG_STRUC.platform > 0) | (ORG_STRUC.numParallelCalcs > 1)
               disp(['JobID=' num2str(POP_STRUC.POPULATION(whichInd).JobID) ]);
            end
            % disp(['check status for step: ' num2str(whichInd)]);               %%%% 检查运行状态
            doneOr = checkStatusC(whichInd);
            % disp(doneOr);
            if doneOr
                % disp(POP_STRUC.POPULATION(whichInd).JobID);
                if POP_STRUC.POPULATION(whichInd).JobID == 0.01   
                   POP_STRUC.POPULATION(whichInd).Step = length([ORG_STRUC.abinitioCode]) + 1;
                elseif POP_STRUC.POPULATION(whichInd).JobID == 0.02   
                   POP_STRUC.POPULATION(whichInd).Step = Step + 1;
                else
                   %%% disp('Reading jobs ...')
                   Error = Reading(ORG_STRUC.abinitioCode(Step),whichInd, indic);
                end

                POP_STRUC.POPULATION(whichInd).JobID = 0;
                if POP_STRUC.POPULATION(whichInd).Error > ORG_STRUC.maxErrors
                    POP_STRUC.POPULATION(whichInd).Done = 1;
                    POP_STRUC.POPULATION(whichInd).ToDo = 0;
                    POP_STRUC.POPULATION(whichInd).Folder=0;
                elseif POP_STRUC.POPULATION(whichInd).Step > length ([ORG_STRUC.abinitioCode])
                    POP_STRUC.POPULATION(whichInd).Done = 1;
                    POP_STRUC.POPULATION(whichInd).ToDo = 0;
                    POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
                    POP_STRUC.POPULATION(whichInd).Number = POP_STRUC.bodyCount;
                    % --- GP-predicted density & uncertainty from gp.csv (best-effort) ---
                    % gp.csv: col1=bodyCount(id), col2=nearest-neighbour index,
                    % col3=residual(distance), col6=density_gp, col7=uncertainty.
                    % Read if this body's row exists; otherwise KEEP inherited value
                    % (keptBest inherits from its parent struct copy, like density/enthalpy).
                    try
                        gpFile = [ORG_STRUC.homePath '/CalcFold1/gp.csv'];
                        % keptBest 晶体是父代的整体 struct 拷贝, 已携带
                        % density_gp/uq/distance。它们不重新走 uspexkit gp,
                        % gp.csv 中无其条目。因此: 字段不存在才初始化为 NaN;
                        % 读不到匹配行时保留已继承的值(与 density/enthalpy
                        % 一样不重算), 而不是覆盖成 NaN。
                        if ~isfield(POP_STRUC.POPULATION(whichInd), 'density_gp')
                            POP_STRUC.POPULATION(whichInd).density_gp = NaN;
                        end
                        if ~isfield(POP_STRUC.POPULATION(whichInd), 'uq')
                            POP_STRUC.POPULATION(whichInd).uq = NaN;
                        end
                        if ~isfield(POP_STRUC.POPULATION(whichInd), 'distance')
                            POP_STRUC.POPULATION(whichInd).distance = NaN;
                        end
                        if exist(gpFile, 'file')
                            gpData = csvread(gpFile, 1, 0);   % skip header row
                            if ~isempty(gpData)
                                bc = POP_STRUC.POPULATION(whichInd).Number;
                                rowIdx = find(gpData(:,1) == bc, 1, 'last');
                                if ~isempty(rowIdx) && size(gpData,2) >= 7
                                    POP_STRUC.POPULATION(whichInd).distance   = gpData(rowIdx,3);
                                    POP_STRUC.POPULATION(whichInd).density_gp = gpData(rowIdx,6);
                                    POP_STRUC.POPULATION(whichInd).uq         = gpData(rowIdx,7);
                                end
                            end
                        end
                    catch
                        % never break the read-back flow on a malformed/missing gp.csv;
                        % 保留已继承的字段值, 不覆盖
                    end
                    % --- end GP density read ---
                    LATTICE     = POP_STRUC.POPULATION(whichInd).LATTICE;
                    COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
                    numIons =     POP_STRUC.POPULATION(whichInd).numIons;
                    MtypeLIST =   POP_STRUC.POPULATION(whichInd).MtypeLIST;
                    numMols = sum(POP_STRUC.POPULATION(whichInd).numMols);
                    atomType = ORG_STRUC.atomType;
                    coordinates = zeros(sum(numMols),3);

                    for i = 1 : sum(numMols)
                        coordinates(i,:)=POP_STRUC.POPULATION(whichInd).MOLECULES(i).MOLCENTER/LATTICE;
                    end
                    
                    %%%% disp('makeMatrices ' ) %%%%%%
                    % --- 快速指纹计算（Cython 加速）---
                    % 第一处：分子质心指纹
                    [order_mol, ~, ~, ~, ok1] = fast_fingerprint(LATTICE, coordinates, numMols, 1);
                    if ok1
                        order = order_mol;
                    else
                        % 回退到 Octave 原版
                        [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, coordinates, numMols, 1);
                        [order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numMols);
                    end
                   
                    for i = 1: sum(numMols)
                        POP_STRUC.POPULATION(whichInd).MOLECULES(i).order = order(i);
                    end
                    
                    % --- 第二处：全原子指纹（带分子内距离过滤）---
                    % 计算 Intra_map 并保存为 .npy 供 Python 使用
                    Intra_map = Intra_MOL_dist(MtypeLIST, numIons, ORG_STRUC.STDMOL);
                    [order, FINGERPRINT, atom_fing, V, ok2, soap_fp] = fast_fingerprint_fullatom(...
                        LATTICE, COORDINATES, numIons, atomType, Intra_map);
                    if ~ok2
                        % 回退到 Octave 原版
                        [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
                        tmp = dist_matrix(1:sum(numIons), :);
                        tmp(find(Intra_map==0)) = 0;
                        dist_matrix(1:sum(numIons),:)  = tmp;
                        [order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
                        soap_fp = [];
                    end
                    %%%%%% disp(order ) %%%%%%
                    POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
                    POP_STRUC.POPULATION(whichInd).struc_entr  = structureQuasiEntropy(whichInd, atom_fing);
                    POP_STRUC.POPULATION(whichInd).S_order     = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);
                    POP_STRUC.POPULATION(whichInd).order       = order;
                    if ~isempty(soap_fp)
                        POP_STRUC.POPULATION(whichInd).soap_fp = soap_fp;
                    end

                    disp('Relaxation is done.')     %%%%%
                    disp(' ')
                    
                    POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
                    WriteIndividualOutput_310(whichInd);
                    POP_STRUC.POPULATION(whichInd).Folder=0;
                end
                safesave ('Current_POP.mat', POP_STRUC)
            else
                disp(['Relaxation is not done for step: ' num2str(whichInd)]); 
            end
        end
    end
end
