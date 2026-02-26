function[cladedNodes] = cladeIt(results);


% snapTime = results.snapTime;
nodeIDs = results.nodeIDs;
cellGroups = results.cellGroups;

cladeCounts = max(cellGroups,[],2);

num = length(nodeIDs);
steps = size(cellGroups,1);
layers = 1:steps;

cladeIDmat = cellGroups + repmat(layers' * num*2,[1 num]);
maxG = max(cellGroups(:));
cladeIDs = unique(cladeIDmat(:));
newGroupIDs = 1:length(cladeIDs);
lookupID(cladeIDs) = newGroupIDs;
cladeIDmat = lookupID(cladeIDmat);

cladeIDs = newGroupIDs;
cladeNum = length(cladeIDs);

clear cladedNodes 
for i = 1:cladeNum
    [y x] = find(cladeIDmat==i); 
    cladedNodes{i} = nodeIDs(x);
end



