function [coor, lat] = Read_GULP_Structure()
fp = fopen('optimized.structure');
running = 1;
count = 1;
dolat = 0;
docoor = 0;

while running
    a = fgetl(fp);
    if findstr(a, 'cell')
        dolat = 1;
    elseif ~(isempty(findstr(a, 'fractional')) & isempty(findstr(a, 'cartesian')))
        docoor = 1;
    elseif dolat==1
        optlat = str2num(a);
        if length(optlat)>3
            optlat(4:6) = optlat(4:6)*pi/180;
            lat = latConverter(optlat);
        else
            lat = 0; 
        end
        dolat=0;
    elseif  docoor == 1
        if findstr(a, 'core')
            tmp = str2num(a(12:50));  %  40 is wrong, should be 50
            coor(count,:) = tmp;
            count = count + 1;
        else
            running = 0;
            docoor=0;
        end
    end
end
fclose(fp);
coor = coor - floor(coor);
