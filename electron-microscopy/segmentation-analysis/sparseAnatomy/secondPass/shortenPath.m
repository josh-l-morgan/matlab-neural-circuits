function[skel] = shortenPath(surfVox,seedPath,tips)

numSurf = size(surfVox.conMat,1);

pathDist = zeros(numSurf,1);
for t = 1:length(tips)
    tip = tips(t);
    tipDist = seedPath.dist(tip);
    for r = 1:numSurf*2
        if tip<1,break,end
       pathDist(tip) = max(pathDist(tip),tipDist); 
       tip =  seedPath.pred(tip);
    end
end


reps = 20;
[path2maxDist] = passMaxAndDist(surfVox,pathDist,reps);

bases = path2maxDist.vox(tips); %find new base for each tip
maxPred = path2maxDist.pred; % combine pred lists
maxPred(maxPred <1) = seedPath.pred(maxPred <1);

skelTip = zeros(1,numSurf);
for t = 1:length(tips)
    tip = tips(t);
    for r = 1:numSurf*2
        if tip<1,break,end
       skelTip(tip) = tips(t); 
       tip =  maxPred(tip);
    end
end

skel.tips = tips;
skel.prop = skelTip;
skel.pred = maxPred;

%%