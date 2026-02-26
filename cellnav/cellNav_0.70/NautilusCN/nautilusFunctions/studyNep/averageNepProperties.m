function[nearVal rec] = averageNepProperties(edges,edgeLength,nodeVal,useNodes,aveDist,recOn);

if ~exist('recOn','var')
    recOn = 0;
end

%aveDist = 25;
nearVal = nodeVal*0;
for i = 1: length(nodeVal);
    %disp(sprintf('%d of %d',i,length(nodeVal)))
    
    pp = node2nodeDist(edges,edgeLength,i,[],aveDist);
    
    nearNodes = find((pp.dists<=aveDist) & useNodes);
    
    nearVals = nodeVal(nearNodes);
    nearVal(i) = mean(nearVals);
    
    if recOn
        recVals{i} = nearVals;
        recDists{i} = pp.dists(nearNodes);
    end
end

if recOn
    rec.vals = recVals;
    rec.dists = recDists;
else
    rec = [];
end
