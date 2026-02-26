function[pathLengths] = itterativeBranchCompression(pred,dists,tips)


%%Draws path distance and owner based on seed paths. By default, distances
%%are incrimented slightly so they can be climbed towards seed
%%by default, all voxels are used as tips.



numSurf = length(dists);

if ~exist('tips','var')
    tips = [1:numSurf];
end


[pathL] = drawPathDist(seedPath.pred,seedPath.dist);
pathDist = floor(pathL.pathDist);
biggestPath = find(pathL.owned == max(pathL.owned),1);
bigPath = pathL.owner == biggestPath;
maxDist = showMaxProp(moveObj,bigPath*1115)+10;
maxDist = showMaxProp(moveObj,pathDist/10)+10;

for i = 1:10
        
        subReps = 5;
        [path2big] = passMaxAndDist(surfVox,bigPath,subReps);
        
        [pathL] = drawPathDist(path2big.pred,path2big.dist);
        pathDist = floor(pathL.pathDist);
        biggestPath = find(pathL.owned == max(pathL.owned),1);
        bigPath = pathL.owner == biggestPath;
        maxDist = showMaxProp(moveObj,bigPath*1115)+10;
        maxDist = showMaxProp(moveObj,pathDist*10)+10;
        

    
end





%% Sort tips
useTips = setdiff(tips,pred);
[sortDists distIDX] = sort( dists(useTips),'descend');
sortTips = useTips(distIDX);


reps = 15;
pathLN = pathL1;
 maxOwned = showMaxProp(moveObj,pathLN.owned+10);
    image(maxOwned(400:600,200:300))
for i = 1:reps
    
    subReps = 5;
    [path2maxDist] = passMaxAndDist(surfVox,pathLN.pathDist,subReps);
    
    newPred = path2maxDist.pred;
    notChanged = sum(newPred<1)
    newPred(notChanged) =  seedPath.pred(notChanged);
    
    %%for each vox on path from tip, enter distance of longest source tip
    [pathLN] = drawPathDist(newPred,seedPath.dist);
    %showMaxProp(moveObj,seedPath.dist/5)
    
    colormap gray(256)
    %showMaxProp(moveObj,pathLN.pathDist+10);
    maxOwned = showMaxProp2(moveObj,pathLN.owned+10);
    maxDist = showMaxProp2(moveObj,pathLN.pathDist/5);
    
    subplot(2,1,1)
    image(maxOwned(400:600,400:600))
    subplot(2,1,2)
        image(maxDist(400:600,400:600))

    
    pause
end