function[path2seed] = conMatRad2allShortestSurfProp(vox);

%

conMat = double(vox.conMat);

    showVal = 1;

for r = 1:26;
    
    %% Grab voxels
    con2 = conMat;
    con3 = con2(:,r);
    isCon = find(con3>0);
    vNum = size(con2,1);
    isSurf = find(con2(:,r)==0);
    seed = isSurf;
    
    %% find shortest path assuming no loops
    
    overShoot = vNum*20;
    dist2seed = ones(1,vNum)*overShoot;
    dist2seed(seed) = 0;
    pred = zeros(vNum,1);
    segID = zeros(vNum,1);
    segID(seed)  = seed;
    dist1 = vox.conDat.dists';
    countBreak = 0;
    distMat = repmat(dist1,[size(con2,1) 1]);
    
    %%
    
    
    nextCon = isSurf; %initialize next verticies to try
    
    
    if showVal
        minI = min(vox.subs,[],1);
        vSubs = vox.subs-repmat(minI,[size(vox.subs,1) 1])+1;
        iSize = max(vSubs,[],1);
        iInd = sub2ind(iSize,vSubs(:,1),vSubs(:,2),vSubs(:,3));
        I = zeros(iSize);
        allI = I;
        allI(iInd) = 1;
        sumI = sum(allI,3);
        I(iInd) = dist2seed;
        minI = sum(I,3);
        colI(:,:,1) = sumI*30;
        colI(:,:,2) = minI;
        colI(:,:,3) = minI;
        
        image(uint8(colI));
        pause(.1)
    end
    
    oldPred = pred;
    oldDist2seed = dist2seed;
    
    for vRep = 1:size(con2,1)
        
        oldDist = dist2seed(isCon)';
        newDist = dist2seed(con3(isCon))' + distMat(isCon,r);
        betterInd = find(oldDist>newDist);
        dist2seed(isCon(betterInd)) = newDist(betterInd);
        pred(isCon(betterInd)) = con3(isCon(betterInd));
        
        segID(isCon(betterInd)) = 1;
        
        if showVal
            
            if ~mod(vRep-1,1)
                vRep
                clf
                showDist = dist2seed;
                
                showDist(showDist>30) = 0;
                
                I(iInd) = showDist;
                minI = sum(I,3);
                colI(:,:,1) = sumI*30;
                colI(:,:,2) = minI*3;
                colI(:,:,3) = minI;
                
                image(uint8(colI));
                %
                %            vRep
                %             %showMaxProp(vox.subs,mod(segID,256)+10);
                %             clf
                %             showAllPred3D(vox.subs,pred);
                %             hold on
                %             %showPred3D(vox.subs,pred);
                pause
            end
        end
        
        if sum(dist2seed) >= sum(oldDist2seed)
            break
        end
        oldPred = pred;
        oldDist2seed = dist2seed;
    end
    
    dist2seed(dist2seed==overShoot) = 0;
    
    d2sMat(:,r) = dist2seed;
    
end

%% analyze d2sMat


shiftMat = ones(3,3,3);
shiftMat(2,2,2) = 0;
shiftInd = find(shiftMat);
[y x z] = ind2sub(size(shiftMat),shiftInd);
vec = [y x z] -2;
op = y*0;
for i = 1:length(y)
    op(i) = find(sum(abs(vec - repmat(vec(i,:)*-1,[size(vec,1) 1])),2) == 0);
end
ori = [[1:13]' op(1:13)];

repel = 1./(d2sMat+1);

for o = 1:26
   oriMat(:,o) = repel(:,o) - repel(:,op(i));
   
end



conD2sMat = d2sMat * 0;
conD2sMat(con2>0) = d2sMat(con2(con2>0));
gradMat = d2sMat-conD2sMat;



%%
path2seed.source = vox.name;
path2seed.seed = seed;
path2seed.pred = pred;
path2seed.dist = dist2seed;
path2seed.segID = segID;
path2seed.d2sMat = d2sMat;
path2seed.repel = repel;
path2seed.oriMat = oriMat; % -1 to 1 13 orientation gradient




