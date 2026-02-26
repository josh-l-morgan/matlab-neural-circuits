function[nep] = sortNodesByLength(nep)

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


if ~(isfield(nep,'pred')) | ~(isfield(nep,'seeds'))
    nep = edges2predFix(nep);
end
seeds = nep.seeds;
pred = nep.pred;

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

%% spread maximum tip dist

tipLength = dists(tips);
[tipLength idx] = sort(tipLength,'descend');
tips = tips(idx);

tipID = nodes*0;
tipID(tips) = 1:length(tips);

for t = 1:length(tips)
    o = tips(t);
    while 1
        n = pred(o);
        if tipID(n)
            break
        else
            tipID(n) = t;
        end
        o = n;
    end
end

nep.tipID = tipID;

if 0% show dist
    colMap = jet(200);
    distInd = round(tipID * 300/max(tipID(:)));
    distInd(distInd<1) = 1;
    distInd(distInd>200) = 200;
    distCol = colMap(distInd,:);
    sP = scatter3(pos(:,1),pos(:,2),pos(:,3),'.');
    sP.CData = distCol;
    drawnow
end




try close(tempFig), end














