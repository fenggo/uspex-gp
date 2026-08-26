function compositionStatistic( composition )
global USPEX_STRUC
global   ORG_STRUC
lenOfBlocks = size(composition, 2);
lenOfComp   = size(composition, 1);
statistic = zeros(lenOfComp, lenOfBlocks*2+7);
statistic(:,1:lenOfBlocks) = composition;
for i = 1:size(composition,1)
statistic(i, lenOfBlocks+1:lenOfBlocks*2) = composition(i,:)/sum( composition(i,:) );
end
for nInd = 1: length(USPEX_STRUC.POPULATION);
if isempty(USPEX_STRUC.POPULATION(nInd).numBlocks)
if ORG_STRUC.molecule ==1
USPEX_STRUC.POPULATION(nInd).numBlocks = USPEX_STRUC.POPULATION(nInd).numMols/ORG_STRUC.numIons;
else
USPEX_STRUC.POPULATION(nInd).numBlocks = USPEX_STRUC.POPULATION(nInd).numIons/ORG_STRUC.numIons;
end
end
whichComposition = [];
numBlocks = USPEX_STRUC.POPULATION(nInd).numBlocks;
for i = 1:lenOfComp
if abs(composition(i,:)-numBlocks) < 1e-3
whichComposition = i;
break;
end
end
if isempty(whichComposition)  
lenOfComp        = lenOfComp +1;
whichComposition = lenOfComp;
composition(lenOfComp,:) =   numBlocks;
statistic(lenOfComp,1:lenOfBlocks*2) = [numBlocks, numBlocks/sum( numBlocks )];
end
statistic(whichComposition, lenOfBlocks*2+1) = statistic(whichComposition, lenOfBlocks*2+1)+1;
switch USPEX_STRUC.POPULATION(nInd).howCome
case '  Random  '
statistic(whichComposition, lenOfBlocks*2+2) = statistic(whichComposition, lenOfBlocks*2+2)+1;
case ' Heredity '
statistic(whichComposition, lenOfBlocks*2+3) = statistic(whichComposition, lenOfBlocks*2+3)+1;
case {'TransMutate','softmutate',' LatMutate ','CoorMutate','Permutate',' Rotate '}
statistic(whichComposition, lenOfBlocks*2+4) = statistic(whichComposition, lenOfBlocks*2+4)+1;
case '  Seeds'
statistic(whichComposition, lenOfBlocks*2+5) = statistic(whichComposition, lenOfBlocks*2+5)+1;
case ' COPEX '
statistic(whichComposition, lenOfBlocks*2+6) = statistic(whichComposition, lenOfBlocks*2+6)+1;
otherwise 
%disp(['Come form : ' USPEX_STRUC.POPULATION(nInd).howCome ])
statistic(whichComposition, lenOfBlocks*2+7) = statistic(whichComposition, lenOfBlocks*2+7)+1;
end
end
USPEX_STRUC.statistic       = statistic;
ORG_STRUC.firstGeneSplitAll = composition;  
