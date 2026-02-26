function[groupNodes] = getSplit(cladeRes,groupNum)


if ~exist('groupNum','var')
    groupNum = 2;
end

nodeIDs = cladeRes.nodeIDs;
cellGroups = cladeRes.cellGroups;

grouped = [];
for i = 1:size(cellGroups,1)
    grouped = cellGroups(i,:);
   if max(grouped) == groupNum;
       break
   end
end


for i = 1:max(grouped);
    groupNodes{i} = nodeIDs(grouped == i);
end