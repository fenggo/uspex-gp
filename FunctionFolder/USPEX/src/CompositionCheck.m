function comp_OK = CompositionCheck(numBlocks)
comp_OK = 1;
if exist('Seeds/Anti-compositions','file')
    [nothing, compStr] = unix('cat Seeds/Anti-compositions |grep -v [[:alpha:]] ');
    comp = str2num( compStr );
    [nothing, str1] = unix('cat Seeds/Anti-compositions |grep s| wc -l ');
    [nothing, str2] = unix('cat Seeds/Anti-compositions |grep S| wc -l ');

    if str2num(str1) || str2num(str2)
        if length( find(numBlocks==0) )==1
            comp_OK = 0;
            return;
        end
    end

    [nothing, str1] = unix('cat Seeds/Anti-compositions |grep b| wc -l ');
    [nothing, str2] = unix('cat Seeds/Anti-compositions |grep B| wc -l ');
    if str2num(str1) || str2num(str2)
        if length( find(numBlocks==0) )==2
            comp_OK = 0;
            return;
        end
    end

    [nothing, str1] = unix('cat Seeds/Anti-compositions |grep t| wc -l ');
    [nothing, str2] = unix('cat Seeds/Anti-compositions |grep T| wc -l ');
    if str2num(str1) || str2num(str2)
        if length( find(numBlocks==0) )==3
            comp_OK = 0;
            return;
        end
    end

    for i = 1:size(comp,1)
        if find( comp(i,:)<0 )                
            if abs( abs(comp(i,:)) - numBlocks ) < 1E-3
                comp_OK = 0;
                return;
            end
        else                                  
            if sameComposition(comp(i,:),numBlocks)
                comp_OK = 0;
                return;
            end
        end
    end
end
