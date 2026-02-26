


springRes = 'D:\LGNs1\Analysis\springDat\results\'
saveName = 'randNoSeedsFilt1res_realNoSeedFilt1Wait30run25.mat'

load([springRes saveName])
results = allResults{1};


    %%
    
snapTime = results.snapTime;
nodeIDs = results.nodeIDs;
cellGroups = results.cellGroups;

cladeCounts = max(cellGroups,[],2);

num = length(nodeIDs);
steps = length(snapTime);
layers = 1:steps;

cladeIDmat = cellGroups + repmat(layers' * num*2,[1 num]);
maxG = max(cellGroups(:));
cladeIDs = unique(cladeIDmat(:));
newGroupIDs = 1:length(cladeIDs);
lookupID(cladeIDs) = newGroupIDs;
cladeIDmat = lookupID(cladeIDmat);
    
    
    
    %%
    cellIDs
    
    
%% structures
maxVal = 1;
for i = 1:length(cladeN)
      clade = cladeN{i};
  maxVal = max(maxVal,max(clade));
end
lookupCheckN = zeros(1,maxVal);
lookupCheckN(checkIDs) = 1:length(checkIDs);


for i = 1:length(cladeN)
   clade = cladeN{i};
   checkPos = lookupCheckN(clade);
   checkedPos = checkPos(checkPos>0);
   L(i) = length(checkedPos); 
end

%% find clade differences

%cladeN = cladeResults.cladedNodes;
allDifs = zeros(sum(L),size(checkProp,2)); allCounts = zeros(sum(L),1);
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
   c = c+countSize;
end

sumDifs = sum(abs(allDifs),2);

difRes.allDifs = allDifs; %difference for each node with a check value from clade mean;
difRes.allCounts = allCounts; % number of checkable nodes in clade
difRes.meanClade = meanClade; %mean value for each property in clade
difRes.sumDifs = sumDifs;
difRes.useDifs = sumDifs(allCounts>1,:);



