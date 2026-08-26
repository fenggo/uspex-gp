function [goodPop, POP, nsym] = PyXtalRandom_Init_310(goodPop, numMols0, nsym)
% PyXtalRandom_Init_310 — PyXtal 替代 Random_Init_310 的初始随机结构生成
%
% 与 Random_Init_310 的区别:
%   - 用 PyXtal random_crystal 替代 symope_310 (random_cell_mol)
%   - 分子放置由 Octave 直接完成，不调用 GetOrientation
%   - 原子顺序由 MOL 文件保证

global ORG_STRUC

badSymmetry = 0;
newSym = 1;
CenterminDist = ORG_STRUC.CenterminDistMatrice;
numMols = numMols0;
running = 1;

while running
    if badSymmetry > 5
        badSymmetry = 0;
        newSym = 1;
    end
    badSymmetry = badSymmetry + 1;

    if newSym
        tmp = find(ORG_STRUC.nsym > 0);
        nsym = tmp(ceil(rand * length(tmp)));
        newSym = 0;
    end

    if ORG_STRUC.constLattice
        targetVol = abs(det(ORG_STRUC.lattice));
    else
        targetVol = ORG_STRUC.latVolume * sum(numMols) / sum(ORG_STRUC.numMols);
    end

    if nsym > 1
        splitInto = ORG_STRUC.splitInto;
        if sum(splitInto) > 3
            numMols = numMols0 / prod(splitInto);
            for i = 1:3
                if size(targetVol) == [1 1]
                    targetVol = targetVol / splitInto(i);
                elseif size(targetVol) == [3 3]
                    targetVol(i, :) = targetVol(i, :) / splitInto(i);
                else
                    targetVol(i) = targetVol(i) / splitInto(i);
                end
            end
        end

        % === 调用 PyXtal ===
        homePath = ORG_STRUC.homePath;
        molStr = '';
        for i = 1:length(ORG_STRUC.numMols)
            if i > 1, molStr = [molStr ',']; end
            molStr = [molStr fullfile(homePath, ['MOL_' num2str(i)])];
        end

        numMolsStr = '';
        for i = 1:length(numMols)
            if i > 1, numMolsStr = [numMolsStr ',']; end
            numMolsStr = [numMolsStr num2str(numMols(i))];
        end

        seed = floor(rand * 2^31);
        pythonBin = '/home/feng/.local/anaconda/bin/python3';
        pyScript = fullfile(homePath, 'FunctionFolder', 'USPEX', '310', 'pyxtal_random.py');
        outDir = fullfile(homePath, 'CalcFoldTemp');
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        cmd = sprintf('%s %s --mols %s --numMols %s --volume %.1f --spg %d --seed %d --outdir %s', ...
            pythonBin, pyScript, molStr, numMolsStr, targetVol, nsym, seed, outDir);
        [status, ~] = system(cmd);

        if status ~= 0
            ORG_STRUC.nsym(nsym) = 0;
            newSym = 1;
            continue;
        end

        dataFile = fullfile(outDir, 'pyxtal_data.txt');
        if ~exist(dataFile, 'file')
            ORG_STRUC.nsym(nsym) = 0;
            newSym = 1;
            continue;
        end

        fid = fopen(dataFile, 'r');
        if fid == -1, continue; end

        header = fscanf(fid, '%d %d', 2);
        nMols = header(1);
        nTypes = header(2);

        lattice = zeros(3, 3);
        for i = 1:3
            lattice(i, :) = fscanf(fid, '%f %f %f', 3);
        end

        molCenters = zeros(nMols, 3);
        molTypes = zeros(nMols, 1);
        for i = 1:nMols
            vals = fscanf(fid, '%f %f %f %f', 4);
            molTypes(i) = vals(1) + 1;
            molCenters(i, :) = vals(2:4);
        end
        fclose(fid);

        % 距离检查
        if distanceCheck(molCenters, lattice, numMols, CenterminDist - 0.2)
            POP.LATTICE = lattice;
            running = 0;
            goodbad = 1;
        else
            goodbad = 0;
        end
    end
end

if goodbad
    [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
    for item = 1:20
        % 放置分子
        lat_6 = latConverter(lattice);
        Molecules = struct('MOLCOORS', {}, 'ZMATRIX', {}, 'ID', {}, 'MOLCENTER', {});
        MtypeLIST2 = zeros(1, nMols);

        for i = 1:nMols
            molType = molTypes(i);
            coords = ORG_STRUC.STDMOL(molType).molecule;
            coords = bsxfun(@minus, coords, mean(coords));

            % 随机旋转
            coords = Rotate_rigid_body([0 0 0], [1 0 0], coords, (rand - 0.5) * 2 * pi);
            coords = Rotate_rigid_body([0 0 0], [0 1 0], coords, (rand - 0.5) * 2 * pi);
            coords = Rotate_rigid_body([0 0 0], [0 0 1], coords, (rand - 0.5) * 2 * pi);

            center_cart = Frac2Cart(molCenters(i, :), lat_6);
            coords = bsxfun(@plus, coords, center_cart);

            Molecules(i).MOLCOORS = coords;
            Molecules(i).MOLCENTER = center_cart;
            format = ORG_STRUC.STDMOL(molType).format;
            Molecules(i).ZMATRIX = NEW_coord2Zmatrix(coords, format);
            MtypeLIST2(i) = molType;
        end

        goodBad = newMolCheck(Molecules, lattice, MtypeLIST2, ORG_STRUC.minDistMatrice);
        if goodBad
            if sum(splitInto) > 3
                [Molecules, numMols, lattice] = SuperMol(Molecules, numMols, lattice, splitInto);
                [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
            end
            [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
            POP.MOLECULES  = Molecules;
            POP.numMols    = numMols;
            POP.MtypeLIST  = MtypeLIST;
            POP.typesAList = typesAList;
            POP.numIons    = numIons;
            POP.LATTICE    = lattice;
            POP.howCome    = '  PyXtal  ';
            goodPop = goodPop + 1;
            disp(['Molecular Crystal ' num2str(goodPop) ...
                  ' built with PyXtal (SPG ' num2str(nsym) ...
                  ' ' spaceGroups(nsym) ')']);
            break;
        end
    end
end
end