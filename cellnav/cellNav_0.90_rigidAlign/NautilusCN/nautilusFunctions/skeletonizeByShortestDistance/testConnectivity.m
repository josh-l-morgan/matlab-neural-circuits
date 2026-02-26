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
    
    nextCon = seed%initialize next verticies to try
    segID(seed) = seed;

    for vRep = 1:size(con2,1)
        vRep
        
        %disp(sprintf('%d ',vRep))
        lastCon = nextCon;
        conTo = con2(lastCon,:); % find who is connected to current vertex
        uniqueIDs = unique(conTo(conTo>0)); %  add connected to list for next start verticies.
        con2(lastCon,:) = 0; % clear connections once vertex is run once

        segID(uniqueIDs) = seed;
        
         showMaxProp(moveObj,mod(segID,256)+10);
        
         nextCon = uniqueIDs;
        
        
        if isempty(nextCon)
            'breaking'
            pause
            break
        end
    end
    
    if sum(segID==0) %% find nearest
       
        countBreak = countBreak + 1
        sourceSub = vox.subs(firstSeed,:);
        emptyVox = find(segID == 0);
        seed = emptyVox(1);
%         
%         targetSub = vox.subs(emptyVox,:);
%         dist2target = sqrt((targetSub(:,1)-sourceSub(1)).^2 + (targetSub(:,2)-sourceSub(2)).^2 ...
%             + (targetSub(:,3)-sourceSub(3)).^2);
%         nearest = find(dist2target == min(dist2target),1);
%         seed = emptyVox(nearest);
%         disp(sprintf('seg ID of new seed is %d',segID(seed)))
        segID(seed) = seed;
        dist2seed(seed) = 0;
        
    end
    
end


