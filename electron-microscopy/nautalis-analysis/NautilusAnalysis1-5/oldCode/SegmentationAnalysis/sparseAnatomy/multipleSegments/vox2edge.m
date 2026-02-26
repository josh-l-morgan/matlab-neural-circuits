%function[pred dist2seed con1] = vox2edge(obj,seed);;

%% study

if ~exist('seed','var')
    seed = 1;
end

vNum = size(obj,1);


%% turn voxel list into connectivity list
con1 = zeros(vNum,26);
dist1 = con1;
conThresh = 1.9;

for i = 1:vNum;
    disp(sprintf('%d of %d',i,vNum))
    tempCon = zeros(1,26);
    tempDist = zeros(1,26);
    dists = sqrt(sum((obj-repmat(obj(i,:),[vNum 1])).^2,2));
    conTo = find((dists >0) & (dists<conThresh));
    tempCon(1:length(conTo)) = conTo;
    con1(i,:) = tempCon;
    tempDist(1:length(conTo)) = dists(conTo);
    dist1(i,:) = tempDist;
end


%% find shortest path assuming no loops
con2 = con1;
dist2seed = ones(1,vNum)*vNum;
dist2seed(seed) = 0;
pred = zeros(vNum,1);

nextCon = seed; %initialize next verticies to try

for vRep = 1:vNum
    lastCon = nextCon;
    nextCon = [];
    for i = 1:length(lastCon) % run all next verticies
        cV = lastCon(i); % grab current vertex
        conTo = con2(cV,:); % find who is connected to current vertex
        con2(cV,:) = 0; % clear connections once vertex is run once
        conTo = conTo(conTo>0); % grab actual connections
        distTo = dist1(cV,conTo>0); %
        betterPath = ( dist2seed(conTo) > (dist2seed(cV)+distTo)); %find connections to which better path has been found
        dist2seed(conTo(betterPath)) = (dist2seed(cV)+1); % Change distances to improve
        pred(conTo(betterPath)) = cV; % mark improved position by predecessor.
        nextCon = cat(2,nextCon,conTo); %  add connected to list for next start verticies.
    end
    
    if ~sum(con2(:,1))
        break
    end
end

