function[loopID] = checkNepForLoops(nep)

% Assign IDs to nodes according to participation in path from furthest tips
% to seed node
%
% global globES
% nep = globES.data.currentNep;
edges = nep.edges;
nodes = nep.nodes;
seed = nep.seedNode;
pos = nep.pos;




if 0
    tempFig = figure;
    scatter(pos(:,1),pos(:,2),'.')
    pause(1)
end


[pred group seeds] = edges2predFix(edges,seed);
%pred = edges2pred(edges,seed);
tips = setdiff(nodes,pred);


%% find distance to seed

predEdge = [[1:length(pred)]' pred];
predLength = sqrt((pos(predEdge(:,1),1) - pos(predEdge(:,2),1)).^2 ...
    + (pos(predEdge(:,1),2) - pos(predEdge(:,2),2)).^2 ...
    + (pos(predEdge(:,1),3) - pos(predEdge(:,2),3)).^2);


dists = nodes * 0  + sum(predLength);

for s = 1:length(seeds)
    newNode = seeds(s);
    dists(newNode) = 0;
    oldDistDif = sum(predLength);
    while 1
        newDists = dists(pred) + predLength';
        dists = min(dists,newDists);
        
        newDistDif = sum(abs(newDists-dists));
        if newDistDif == oldDistDif
            break
        end
        oldDistDif = newDistDif;
        
        if 0% show dist
            colMap = jet(200);
            distInd = round(dists * 200/max(dists(:)));
            distInd(distInd<1) = 1;
            distInd(distInd>200) = 200;
            distCol = colMap(distInd,:);
            sP = scatter3(pos(:,1),pos(:,2),pos(:,3),'.');
            sP.CData = distCol;
            drawnow
            
        end
        
        
    end
    
    if oldDistDif
        
        disp('Poorly conditioned skeleton')
        
    end
end