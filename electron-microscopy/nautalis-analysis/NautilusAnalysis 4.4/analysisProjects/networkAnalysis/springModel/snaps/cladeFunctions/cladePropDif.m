function[difRes] = cladePropDif(cladeN,checkIDs,checkProp)

%%CladedNodes is a cell array where each cell is a clade and the contents
%%are the nodes included in the clade


%% structures
maxVal = 1;
for i = 1:length(cladeN)
      clade = cladeN{i};
  maxVal = max(maxVal,max(clade));
end
lookupCheckN = zeros(1,maxVal);
lookupCheckN(checkIDs) = 1:length(checkIDs); %value at position is check index


for i = 1:length(cladeN)
   clade = cladeN{i};
   checkPos = lookupCheckN(clade);
   checkedPos = checkPos(checkPos>0);
   L(i) = length(checkedPos); 
end

checkSlots = [0 cumsum(L)];
allCheckedVals = zeros(max(checkSlots),size(checkProp,2));
for i = 1:length(cladeN)
   clade = cladeN{i};
   checkPos = lookupCheckN(clade);
   checkedPos = checkPos(checkPos>0);
   allCheckedVals(checkSlots(i)+1:checkSlots(i+1),:) = checkProp(checkedPos,:);
end
meanVal = mean(allCheckedVals,1);

%% find clade differences

%cladeN = cladeResults.cladedNodes;
allDifs = zeros(sum(L),size(checkProp,2)); allCounts = zeros(sum(L),1);
allNormDifs = allDifs; allNormCor = allDifs;
Lsum = [1 cumsum(L)];
c = 1;
for i = 1:length(cladeN)
    
   clade = cladeN{i};
   checkPos = lookupCheckN(clade);
   checkedPos = checkPos(checkPos>0);
   cladeProps = checkProp(checkedPos,:);
   difProp = cladeProps - repmat(mean(cladeProps,1),[size(cladeProps,1),1]);
    countSize = length(checkedPos);
  
   allDifs(c+1: c+countSize,:) = difProp;
   allCounts(c+1: c+countSize) = repmat(countSize,[countSize 1]);
   meanClade(i,:) = mean(cladeProps,1);
   
   
   
   normProp = cladeProps - repmat(meanVal,[size(cladeProps,1) 1]);
   meanNorm = mean(normProp);
   normDifs = normProp - repmat(meanNorm,[size(normProp,1),1]);
   normCor = normProp .* repmat(meanNorm,[size(normProp,1),1]);
   allNormDifs(c+1: c+countSize,:) = normDifs;
   allNormCor(c+1: c+countSize,:) = normCor;
      
   c = c+countSize;
end

sumDifs = sum(abs(allDifs),2);

difRes.allDifs = allDifs; %difference for each node with a check value from clade mean;
difRes.allCounts = allCounts; % number of checkable nodes in clade
difRes.meanClade = meanClade; %mean value for each property in clade
difRes.sumDifs = sumDifs;
difRes.allNormCor = allNormCor;
difRes.allNormDfi = allNormDifs;
difRes.wcss = sum(allDifs.^2,2);
difRes.useDifs = sumDifs(allCounts>1,:);
difRes.useWCSS = difRes.wcss(allCounts>1,:);
difRes.useNC = allNormCor(allCounts>1,:);









