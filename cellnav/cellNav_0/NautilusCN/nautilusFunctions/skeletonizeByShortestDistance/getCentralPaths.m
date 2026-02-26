function[centralPath] = getCentralPaths(surfVox)


seedPath = surfVox.seedPath;

voxList = ones(size(surfVox.subs,1),1);
usedPred = unique(seedPath.pred);
usedPred = usedPred(usedPred>0);
voxList(usedPred) = 0;
tips = find(voxList);

numSurf = size(surfVox.conMat,1);

%% for each vox on path from tip, enter distance of longest source tip


pathDist = zeros(numSurf,1);
for t = 1:length(tips)
    tip = tips(t);
    tipDist = seedPath.dist(tip);
    for r = 1:numSurf*2
        if tip<1,tip;break,end
        pathDist(tip) = max(pathDist(tip),tipDist);
        tip =  seedPath.pred(tip);
    end
end

%showMaxProp(moveObj,pathDist+ones(size(moveObj,1),1)*10);

%% find path to voxel that is path to greatest distnce tip
reps = 30;
[path2maxDist] = passMaxAndDist(surfVox,pathDist,reps);

useOld = path2maxDist.pred<1;
bases = path2maxDist.vox; %find new base for each tip
%bases(useOld) = seedPath.vox(useOld);
maxPred = path2maxDist.pred; % combine pred lists
maxPred(useOld) = seedPath.pred(useOld); %keep preds with no new solutions
minDist = path2maxDist.dist;
minDist(useOld) = seedPath.dist(useOld);

%% for each vox on path from tip, enter distance of longest source tip


newPathDist = zeros(numSurf,1);
for t = 1:length(tips)
    tip = tips(t);
    tipDist = minDist(tip);
    for r = 1:numSurf*2
        if tip<1,tip;break,end
        newPathDist(tip) = max(newPathDist(tip),tipDist);
        tip =  maxPred(tip);
    end
end

%showMaxProp(moveObj,newPathDist+ones(size(moveObj,1),1)*10);

%% Record
centralPath.tips = tips;
centralPath.bases = bases;
centralPath.dists = minDist;
centralPath.pred = maxPred;
centralPath.pathDist = newPathDist;


