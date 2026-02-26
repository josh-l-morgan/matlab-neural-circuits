function[pathD] = setDist2Seed(pathD)
%%Recalculate distance to seed




pred = pathD.pred(:);
predL = pathD.predLength(:);
dists = pred * 0;
currentPred = pred;

for i = 1:length(pred)
   
    doShift = find(currentPred>0);
    if isempty(doShift)
        break
    end    
    dists(doShift) = dists(doShift) + predL(currentPred(doShift));
    currentPred(doShift) = pred(currentPred(doShift));
%     [sortPaths pIDX]  = sort(pathD.lengths,'descend');
%     showTips = pIDX(1:50);
%     showPathL3D3(allVox.subs,propPathI,dists,showTips)
%     pause(.1)

end

pathD.lengths = dists(pathD.tips)';
pathD.voxLengths = dists(pathD.owner);
pathD.dist = dists;

pathD = drawPathDist(pathD.pred,pathD.dist,pathD.predLength,pathD.tips);  %% organize all paths according furthest from seed
































