function[skel] = shortenPath(surfVox)


seedPath = surfVox.seedPath;
tips = surfVox.longTip.tips;
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

%% Create bones
for t = 1:length(tips)
    tip = tips(t); 
    bones(t).tip = tip;
    boneLength = path2maxDist.dist(tip);
    if ~boneLength
        boneLength = pathDist(tip);
    end
    bones(t).length = boneLength;
    nodes = [];
    for r = 1:numSurf*2
        if tip<1,
            bones(t).base = nodes(end);
            bones(t).parent = 0;
            break,
        end
        
        if skelTip(tip)
            if pathDist(skelTip(tip)) > pathDist(tips(t)) %if new path is from  shorter process
                bones(t).base = tip;
                bones(t).parent = skelTip(tip);
                break
            else
                nodes = [nodes tip];
            end
        else
            nodes = [nodes tip];
        end
        skelTip(tip) = tips(t);
        tip =  maxPred(tip);
        
    end
    bones(t).nodes = nodes;
end

skel.tips = tips;
skel.prop = skelTip;
skel.pred = maxPred;
skel.bones = bones;

%%