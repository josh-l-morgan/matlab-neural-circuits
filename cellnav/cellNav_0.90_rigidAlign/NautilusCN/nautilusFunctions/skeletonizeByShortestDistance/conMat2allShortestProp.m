function[path2seed] = conMat2allShortestProp(vox,prop,reps);

%% study
tic
if ~exist('firstSeed','var')
    firstSeed = 1;
end

if ~exist('reps','var')
    reps = size(vox.conMat,2)*2;
end


ifShow = 0;

seed = firstSeed;

%% Grab voxels
con2 = vox.conMat;
vNum = size(con2,1);
vNum2 = vNum+1;
con2(con2==0) = vNum2;

penalty = prop;
penalty(vNum2) = inf;

%% find shortest path assuming no loops

overShoot = vNum*20;
dist2seed = ones(1,vNum2)*0;  %all prop voxels start at zero
dist2surf = penalty;%
dist2surf(vNum2) = inf;
dist2surfPath = dist2surf;
sourceVox = [1:vNum2];
%dist2seed(seed) = 0;
dist2seed(vNum2) = inf;
pred = zeros(vNum,1);
predLength = zeros(vNum,1);
segID = zeros(vNum,1);
segID(seed)  = seed;
dist1 = vox.conDat.dists';
countBreak = 0;

distMat = repmat(dist1,[size(con2,1) 1]);
%%%
oldPred = pred;
oldDist2seed = dist2seed;
checkV = 14;

    nextCon = seed; %initialize next verticies to try

    %% test data
    %{
    tNum = 4;
    dist2seed = zeros(1,tNum+1);
    dist2seed(end) = inf;
    dist2surfPath = 1:tNum+1;
    dist2surfPath(end) = inf;
    con2 = zeros(tNum,26)+tNum+1;
    con2(3,1) = 1;
    distMat = repmat(dist1,[tNum 1]);
    pred = dist2seed*0;
    sourceVox = [1:tNum];
    %}
    
    
    
    
    
    %%
    for vRep = 1:reps
       
        oldDist = dist2seed(1:end-1)';
        oldValue = dist2surfPath(1:end-1)';
        newDist = dist2seed(con2) + distMat;
        newValue = dist2surfPath(con2);
        bestValue = repmat(min( min(newValue,[],2),oldValue),[1 26]); %must be this value to be of interest
        isBestValue = newValue == bestValue; % find directions were value is equal to best possible value
        isBetterValue = newValue < repmat(oldValue,[1 26]);
        %betterMove = isBestValue & isBetterValue;%
        valDist = newDist;
        %%Move if equal and closer or if better
        valDist(~isBestValue) = inf; % find shorter distance to best new value, or find shorter distance to equal value if none or better
        bestNewValDist = repmat(min(valDist,[],2),[1 26]); %must be this value to be of interest
        isBestValDist = (valDist  == bestNewValDist); % closest directions where best value is better
        
        foundBetterVal = repmat(sum(isBetterValue,2),[1 26])>0; %matrix for zeroing out pure distance if better value is found       
        betterDist = (repmat(oldDist,[1 26])>newDist) & (repmat(oldValue,[1 26]) == newValue); %  better distance and equal found val
        
        betterInd =  find((betterDist & ~foundBetterVal) | (isBestValDist & isBetterValue));
        
        [sourceCon xnul] = ind2sub(size(newDist),betterInd);
        dist2seed(sourceCon) = newDist(betterInd);
        pred(sourceCon) = con2(betterInd);
        dist2surfPath(sourceCon) = dist2surfPath(con2(betterInd)); %inheret prop value of pred
        sourceVox(sourceCon) = sourceVox(con2(betterInd));
        predLength(sourceCon) = distMat(betterInd);
        
        
        segID(sourceCon) = seed;
        
        if ifShow
            if ~mod(vRep-1,1)
                %dist2surfPath(checkV)
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
    
%     if sum(segID==0) %% find nearest
%         countBreak = countBreak + 1
%         sourceSub = vox.subs(firstSeed,:);
%         emptyVox = find(segID == 0);
%         targetSub = vox.subs(emptyVox,:);
%         dist2target = sqrt((targetSub(:,1)-sourceSub(1)).^2 + (targetSub(:,2)-sourceSub(2)).^2 ...
%             + (targetSub(:,3)-sourceSub(3)).^2);
%         nearest = find(dist2target == min(dist2target),1);
%         seed = emptyVox(nearest);
%         disp(sprintf('seg ID of new seed is %d',segID(seed)))
%         segID(seed) = seed;
%         dist2seed(seed) = 0;
%         
%         
%     end
%     
% end


dist2seed(dist2seed==overShoot) = 0;

%%
path2seed.source = vox.name;
path2seed.seed = seed;
path2seed.pred = pred(1:vNum);
path2seed.dist = dist2seed(1:vNum);
path2seed.segID = segID(1:vNum);
path2seed.prop = dist2surfPath(1:vNum);
path2seed.vox = sourceVox(1:vNum);
path2seed.predLength = predLength(1:vNum);
