function[path2seed] = conMat2allShortest(vox,firstSeed);

%% study
tic
if ~exist('firstSeed','var')
    firstSeed = 1;
end

ifShow = 1;

seed = firstSeed;

%% Grab voxels
con2 = vox.conMat;
vNum = size(con2,1);
vNum2 = vNum+1;
con2(con2==0) = vNum2;



penalty = vox.path2surf.penalty;
penalty(vNum2) = inf;

%% find shortest path assuming no loops

overShoot = vNum*20;
dist2seed = ones(1,vNum2)*overShoot;
dist2seed(seed) = 0;
dist2seed(vNum2) = inf;
pred = zeros(vNum,1);
segID = zeros(vNum,1);
segID(seed)  = seed;
dist1 = vox.conDat.dists';
countBreak = 0;

distMat = repmat(dist1,[size(con2,1) 1]);
%%%
oldPred = pred;
oldDist2seed = dist2seed;
while sum(segID==0)
    
    nextCon = seed; %initialize next verticies to try

    for vRep = 1:size(con2,1)
       
        oldDist = dist2seed(1:end-1)';
        newDist = dist2seed(con2) + distMat;
        betterPath = repmat(oldDist,[1 26])>(newDist + penalty);
        minPath = min(newDist,[],2);
        isMin = (newDist == repmat(minPath,[1 26]));
        hit = isMin & betterPath;
        
        sum(sum(hit,2)>1);
        
        betterInd = find(hit); %connected voxels with better paths
        [sourceCon xnul] = ind2sub(size(betterPath),betterInd);
        dist2seed(sourceCon) = newDist(betterInd);
        pred(sourceCon) = con2(betterInd);
        
        segID(sourceCon) = seed;
        
           if ifShow     
       if ~mod(vRep-1,1)
           
           vRep
            %showMaxProp(vox.subs,mod(segID,256)+10);
            clf
            %showAllPred3D(vox.subs,pred);
            path.pred = pred;
            path.dist = dist2seed;
            path.seed = seed;
            showLongPred3D(vox.subs,path);
             hold on
            pause(.01)
       end
           end
        
        if sum(oldPred~=pred)==0
            break
        end
        oldPred = pred;
        oldDist2seed = dist2seed;
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
        disp(sprintf('seg ID of new seed is %d',segID(seed)))
        segID(seed) = seed;
        dist2seed(seed) = 0;
        
        
    end
    
end


dist2seed(dist2seed==overShoot) = 0;

%%
path2seed.source = vox.name;
path2seed.seed = seed;
path2seed.pred = pred(1:vNum);
path2seed.dist = dist2seed(1:vNum);
path2seed.segID = segID(1:vNum);




