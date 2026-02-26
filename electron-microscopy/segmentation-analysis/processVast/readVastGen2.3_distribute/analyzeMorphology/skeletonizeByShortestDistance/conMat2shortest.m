%function[pred ] = conMat2shortests(conMat,conDat,seed);

%% study
tic
if ~exist('seed','var')
    seed = 1;
end

vNum = size(conMat,1);

%% find shortest path assuming no loops
con2 = conMat;
dist2seed = ones(1,vNum)*vNum*2;
dist2seed(seed) = 0;
pred = zeros(vNum,1);
dist1 = conDat.dists';
nextCon = seed; %initialize next verticies to try

for vRep = 1:vNum
    
    disp(vRep)
    tic
    lastCon = nextCon;
    conTo = con2(lastCon,:); % find who is connected to current vertex
    %nextCon = unique(conTo(conTo>0)); %  add connected to list for next start verticies.
    
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
    toc
    %%resolve conflicts
    tic
    uniqueIDs = unique(changeID); %consider using for nextCon
    nextCon = uniqueIDs;
    sourceID = uniqueIDs * 0;
    parfor u = 1:length(uniqueIDs)
        cu = find(changeID == uniqueIDs(u));
        smallDists = find(changeDist(cu) == min(changeDist(cu)));
        sourceID(u) = cu(smallDists(ceil(rand*length(smallDists))));
    end
    toc
   
    pred(uniqueIDs) = predList(sourceID);
    dist2seed(uniqueIDs) = changeDist(sourceID);
   
    if isempty(nextCon)
        break
    end
end


toc