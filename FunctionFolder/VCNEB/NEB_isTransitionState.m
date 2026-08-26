function isTS=NEB_isTransitionState()
global ORG_STRUC
global POP_STRUC
numImages = ORG_STRUC.numImages;
deltaE=0;
for i = 1:numImages
H(i)=POP_STRUC.POPULATION(i).Enthalpy;
end
isTS=zeros(1,numImages);  
for i = 2:numImages-1
if (H(i)-H(i-1)>deltaE) & (H(i)-H(i+1)>deltaE)
isTS(i)=1;
end
end
for i = 2:numImages-1
if (H(i)-H(i-1)<-deltaE) & (H(i)-H(i+1)<-deltaE)
isTS(i)=-1;
end
end
seq_TSH=find(isTS==1);
[seq_TS, result] = sort( H(find(isTS== 1)));
for i = 1:length(seq_TS);
isTS(seq_TSH(i))=i; 
end
seq_LMH=find(isTS==-1);
[seq_LM, result] = sort(-H(find(isTS==-1)));
for i = 1:length(seq_LM);
isTS(seq_LMH(i))=-i; 
end
