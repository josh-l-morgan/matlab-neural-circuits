%function[path2seed] = conMat2allShortest(vox,firstSeed);

%% study
tic
if ~exist('firstSeed','var')
    firstSeed = 1;
end

seed = firstSeed;

%% Grab voxels
con2 = vox.conMat;
vNum = size(con2,1);

%% find shortest path assuming no loops
overShoot = vNum*20;
dist2seed = ones(1,vNum)*overShoot;
dist2seed(seed) = 0;
pred = zeros(vNum,1);
segID = zeros(vNum,1);
segID(seed)  = seed;
dist1 = vox.conDat.dists';
countBreak = 0;

%%

while sum(segID==0)
    
    nextCon = seed; %initialize next verticies to try
    
    for vRep = 1:size(con2,1)
        
        %disp(sprintf('%d ',vRep))
        lastCon = nextCon;
        conTo = con2(lastCon,:); % find who is connected to current vertex
        %nextCon = unique(conTo(conTo>0)); %  add connected to list for next start verticies.
        deepVox = sum(conTo(:,1:6)>0,2)>5; % find deep voxels
        conTo(deepVox,:) = 0; % clear connections of deep voxels
        con2(lastCon,:) = 0; % clear connections once vertex is run once
        conTo(conTo==0) = seed;
        oldDist = dist2seed(lastCon)';
        currentDist = dist2seed(conTo);
        potentialDist = repmat(oldDist,[1 26]) + repmat(dist1,[length(lastCon) 1]);
        betterPath = currentDist>potentialDist;
        
        betterInd = find(betterPath);
        [sourceCon xnul] = ind2sub(size(betterPath),betterInd);
        predList = lastCon(sourceCon);
        changeID = conTo(betterInd);
        changeDist = potentialDist(betterInd);
        
        
        uniqueIDs = unique(changeID); %consider using for nextCon
        nextCon = uniqueIDs;
        sourceID = uniqueIDs * 0;
        parfor u = 1:length(uniqueIDs)
            cu = find(changeID == uniqueIDs(u));
            smallDists = find(changeDist(cu) == min(changeDist(cu)));
            sourceID(u) = cu(smallDists(ceil(rand*length(smallDists))));
        end
        
        pred(uniqueIDs) = predList(sourceID);
        dist2seed(uniqueIDs) = changeDist(sourceID);
        segID(uniqueIDs) = seed;
        
        if isempty(nextCon)
            break
        end
    end
    
    if sum(segID==0) %% find nearest
        countBreak = countBreak + 1
        sourceSub = vox.subs(firstSeed,:);
        emptyVox = find(segID == 0);
        targetSub = vox.subs(emptyVox,:);
        dist2target = sqrt((targetSub(:,1)-sourceSub(1)).^2 + (targetSub(:,2)-sourceSub(2)).^2 ...
            + (targetSub(:,3)-sourceSub(3)).^2);
        nearest = find(dist2target == min(dist2target),1);
        seed = emptyVox(nearest);
        segID(seed) = seed;
        dist2seed(seed) = 0;
        
    end
    
    
end


dist2seed(dist2seed==overShoot) = 0;

%%
path2seed.source = vox.name;
path2seed.seed = seed;
path2seed.pred = pred;
path2seed.dist = dist2seed;
path2seed.segID = segID;




