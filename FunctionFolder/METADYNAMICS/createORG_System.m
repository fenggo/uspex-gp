function createORG_System(inputFile)
    global ORG_STRUC
    
    % 获取 Python 脚本路径
    getPy = [ORG_STRUC.USPEXPath, '/FunctionFolder/getInput.py'];
    
    % 读取计算类型参数
    calculationType = python_uspex(getPy, ['-f ' inputFile ' -b calculationType -c 1']);
    if ~isempty(calculationType)
        ORG_STRUC.dimension = str2num(calculationType(1));
        ORG_STRUC.molecule  = str2num(calculationType(2));
        ORG_STRUC.varcomp   = str2num(calculationType(3));
    end
    
    % 读取 pickUpYN 参数
    pickUpYN = python_uspex(getPy, ['-f ' inputFile ' -b pickUpYN -c 1']);
    if isempty(pickUpYN)
        pickUpYN = '0'; 
    end
    ORG_STRUC.pickUpYN = str2num(pickUpYN);
    
    % 读取 pickUpGen 参数
    pickUpGen = python_uspex(getPy, ['-f ' inputFile ' -b pickUpGen -c 1']);
    if isempty(pickUpGen)
        pickUpGen = '1'; 
    end
    ORG_STRUC.pickUpGen = str2num(pickUpGen);
    
    % 如果 pickUpYN 为 0，强制 pickUpGen 为 1
    if ORG_STRUC.pickUpYN == 0
        ORG_STRUC.pickUpGen = 1; 
    end
    
    % 读取 pickUpFolder 参数
    pickUpFolder = python_uspex(getPy, ['-f ' inputFile ' -b pickUpFolder -c 1']);
    if isempty(pickUpFolder)
        pickUpFolder = '1'; 
    end
    ORG_STRUC.pickUpFolder = str2num(pickUpFolder);
    
    % 打开 POSCAR 文件读取结构信息
    [fid, message] = fopen('POSCAR_1');
    tmp = fgetl(fid);
    scale_factor = fgetl(fid);
    lat = fscanf(fid, '%g', [3, 3]);  % 晶格矢量
    lat = lat';
    tmp = fgetl(fid);
    atomType = fgetl(fid);
    numIons = fgetl(fid);
    
    % 判断是否为 Direct 坐标格式
    if strcmp(numIons, 'Direct')
        ORG_STRUC.numIons = str2num(atomType);
        atomType = python_uspex(getPy, ['-f ' inputFile ' -b atomType -e EndAtomType']);
        atomType(end) = [];
    else    
        ORG_STRUC.numIons = str2num(numIons);
        tmp = fgetl(fid);
    end
    
    % 解析原子类型
    ORG_STRUC.atomType = zeros(1, size(ORG_STRUC.numIons, 2));
    c1 = findstr(atomType, ' ');
    c = sort(str2num(['0 ' num2str(c1)]));
    c(end+1) = length(atomType) + 1;
    ind1 = 1;
    
    for i = 2 : length(c)
        if c(i-1)+1 > c(i)-1
            continue;
        end
        tmp = atomType(c(i-1)+1 : c(i)-1);
        
        if ~isempty(str2num(tmp))
            ORG_STRUC.atomType(ind1) = str2num(tmp);
        else
            for j = 1 : 105
                if strcmp(lower(tmp), lower(elementFullName(j))) | strcmp(lower(tmp), lower(megaDoof(j)))
                    ORG_STRUC.atomType(ind1) = j;
                    break;
                end
            end
        end
        ind1 = ind1 + 1;
    end
    
    % 读取原子坐标
    sss = fscanf(fid, '%g', [3, sum(ORG_STRUC.numIons)]);
    ss = sss';
    coordinates = ss(:, 1:3);
    fclose(fid);
    
    % 计算基矢量
    base(1, :) = lat(1, :) / norm(lat(1, :));
    base(3, :) = cross(lat(1, :), lat(2, :));
    base(3, :) = base(3, :) / norm(base(3, :));
    base(2, :) = cross(base(3, :), base(1, :));
    
    % 存储晶格和坐标
    ORG_STRUC.lattice = lat / base;
    ORG_STRUC.coordinates = coordinates - floor(coordinates);
    
    % 初始化成键矩阵
    ORG_STRUC.goodBonds = zeros(length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        for j = i : length(ORG_STRUC.atomType)
            ORG_STRUC.goodBonds(i, j) = 0.15;
            ORG_STRUC.goodBonds(j, i) = 0.15;
        end
    end
    
    % 读取价电子信息
    valences = python_uspex(getPy, ['-f ' inputFile ' -b valences -e endValences']);
    ORG_STRUC.valences = str2num(valences);
    
    if isempty(valences)
        ORG_STRUC.valences = zeros(1, length(ORG_STRUC.atomType));
        for i = 1 : length(ORG_STRUC.atomType)
            ORG_STRUC.valences(i) = str2num(valence(ORG_STRUC.atomType(i)));
        end
    else
        ORG_STRUC.valences = str2num(valences);
    end
    
    % 读取价电子数
    NvalElectrons = python_uspex(getPy, ['-f ' inputFile ' -b valenceElectr -e endValenceElectr']);
    
    if isempty(NvalElectrons)
        ORG_STRUC.NvalElectrons = zeros(1, length(ORG_STRUC.atomType));
        for i = 1 : length(ORG_STRUC.atomType)
            ORG_STRUC.NvalElectrons(i) = str2num(valenceElectronsNumber(ORG_STRUC.atomType(i)));
            if ORG_STRUC.NvalElectrons(i) == 0
                ORG_STRUC.NvalElectrons(i) = ORG_STRUC.valences(i);
            end
        end
    else
        ORG_STRUC.NvalElectrons = str2num(NvalElectrons);
    end
    
    % 读取离子距离参数
    hardCore = python_uspex(getPy, ['-f ' inputFile ' -b IonDistances -e EndDistances']);
    
    if isempty(hardCore)
        ORG_STRUC.minDistMatrice = zeros(1, length(ORG_STRUC.atomType));
        for i = 1 : length(ORG_STRUC.atomType)
            s = covalentRadius(ORG_STRUC.atomType(i));
            ORG_STRUC.hardCore(i) = str2num(s) / 2;
        end
    else
        ORG_STRUC.hardCore = str2num(hardCore);
    end
    
    % 构建最小距离矩阵
    ORG_STRUC.minDistMatrice = zeros(length(ORG_STRUC.atomType), length(ORG_STRUC.atomType));
    
    if size(ORG_STRUC.hardCore, 1) == size(ORG_STRUC.hardCore, 2)
        for i = 1 : length(ORG_STRUC.atomType)
            for j = i : length(ORG_STRUC.atomType)
                ORG_STRUC.minDistMatrice(i, j) = ORG_STRUC.hardCore(i, j);
                ORG_STRUC.minDistMatrice(j, i) = ORG_STRUC.hardCore(i, j);
            end
        end
    else
        ORG_STRUC.minDistMatrice = zeros(length(ORG_STRUC.hardCore));
        for hardInd_1 = 1 : length(ORG_STRUC.hardCore)
            for hardInd_2 = hardInd_1 : length(ORG_STRUC.hardCore)
                ORG_STRUC.minDistMatrice(hardInd_1, hardInd_2) = ...
                    ORG_STRUC.hardCore(hardInd_1) + ORG_STRUC.hardCore(hardInd_2);
                ORG_STRUC.minDistMatrice(hardInd_2, hardInd_1) = ...
                    ORG_STRUC.hardCore(hardInd_1) + ORG_STRUC.hardCore(hardInd_2);
            end
        end
    end
    
    % 设置最优基础结构
    ORG_STRUC.bestBasicStructure = ORG_STRUC.coordinates;
end
