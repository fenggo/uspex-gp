function target = Read_GULP(flag, ID)
if     flag == -1
    [nothing, results] = unix('grep Finished output');
    if isempty(results)
    USPEXmessage(1151,'',0);
    target = 0;
    else
    target = 1;
    end
elseif flag == 0
    target = 1;
    [nothing, results] = unix(['./getStuff output Energy 5 |tail -1']);
    if isempty(results) | isempty(str2num(results))
        USPEXmessage(1152,'',0);
        target = 0;
    else
        [nothing, results3] = unix(['grep dump optimized.structure']);
        if isempty(results3)
        target = 0;
        end
    end
elseif flag ==  1 
    [nothing, Energy_Str] = unix(['./getStuff output Energy 5 |tail -1']);
    target = str2num(Energy_Str);
elseif flag ==  2 
    target = callAWK('GULP_pres.awk','output');
elseif flag ==  3 
    d_s = zeros(1,6);
    d_s_ion = zeros(1,6);
    [fileHandle, msg] = fopen('output');

    while 1
        tline = fgetl(fileHandle);
        if ~ischar(tline)
        break
        end
        if ~isempty(strfind(tline, 'Static dielectric constant'))
            tline = fgetl(fileHandle);
            tline = fgetl(fileHandle);
            tline = fgetl(fileHandle);
            tline = fgetl(fileHandle);
            tline1 = fgetl(fileHandle);
            e_tensor = sscanf(tline1,'%*s %g %g %g');
            d_s(1) = e_tensor(1);
            d_s(4) = e_tensor(2);
            d_s(5) = e_tensor(3);
            tline2 = fgetl(fileHandle);
            e_tensor = sscanf(tline2,'%*s %g %g %g');
            d_s(2) = e_tensor(2);
            d_s(6) = e_tensor(3);
            tline3 = fgetl(fileHandle);
            e_tensor = sscanf(tline3,'%*s %g %g %g');
            d_s(3) = e_tensor(3);
        end
    end
    fclose(fileHandle);
    target = d_s;
elseif flag ==  6 
    elas_orig = callAWK('GULP_elas.awk', 'output');   
    try 	
        elasticMatrix = elas_orig([1,2,3,4,5,6],:);
        catch
        %disp('Warning : Cannot extract Elastic Constant Matrix from GULP output file! ');
        USPEXmessage(1155, '', 0);            
        elasticMatrix=[];
    end
    target = elasticMatrix;	
elseif  flag ==  7 
    try
        [nothing, numStr] = unix(['grep "Total number atoms/shells" output |cut -d"=" -f2']);
        numIons= sum( str2num(numStr) );
        force_orig = callAWK('GULP_force.awk', 'output', ['num=', num2str(numIons)]);
        if isempty(force_orig)
           USPEXmessage(1156,'',0);
           target = [];
        else
           target = force_orig;
        end
        catch
        target = [];
    end
end
