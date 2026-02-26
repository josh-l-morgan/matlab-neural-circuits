%%Make predictions of how many synapses each pair of cells should form
%%based of the distribution of their skeletons. -
%-Skeletons are filtered



%% Plot axon to tcr
clear all
%MPN = GetMyDir;
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

seedList = [108 201 109 903 907]
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
preList = useList.preList;
colormap gray(256)


%%
voxLength = 0.2; %size of the voxels to be generated for measuring volumes
voxVol = voxLength^3; 
anchorScale =  (  [obI.em.res].*[4 4 1])./ (1000*voxLength )   


%voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3)* 4 ];
voxelScale = obI.em.dsRes;

butSize.voxLength = voxLength;
butSize.voxVol = voxVol;
butSize.anchorScale = anchorScale;
butSize.voxelScale = voxelScale;
%% Get cells

axList = preList;
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

butSize.edges = edges;

rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
rawSynAnchors(rawSynAnchors<1) = 1;
synAnchors = scaleSubs(double(rawSynAnchors),anchorScale);

maxRad = 20;
clear balls
for r = 1:maxRad;
    ball = ones(r*2,r*2,r*2);
    ballInd = find(ball);
    [y x z] = ind2sub(size(ball),ballInd);
    
    dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+(z-mean(z)).^2);
    ball(dists>r) = 0;
    ball = ball/sum(ball(:));
    image(squeeze(sum(ball,2))*100000);
    balls{r} = ball;
end

%%  Get Skeletons
colormap jet(256)
checkRad = maxRad*2;
butSize.axList = axList;
for a = 1:length(axList)
    disp(sprintf('running %d of %d',a,length(axList)))
    %subs = scaleSubs(getCellSubs(obI,dsObj,axList(a)),voxelScale);
        subs = getCellSubs(obI,dsObj,axList(a));

    isAx = find((edges(:,2) == axList(a)) & sum(synAnchors,2));
    synDat.isAx = isAx;
    synDat.edges = edges(isAx,:);
    synDat.anchors = synAnchors(isAx,:);
    synDat.rawAnchors = rawAnchors(isAx,:);
    butSize.synDat(a) = synDat;
    
    scatter(subs(:,1),subs(:,2),'.')
    hold on
    scatter(synAnchors(isAx,1),synAnchors(isAx,2),'r','.')
    hold off
    sc = 0;
    clear butProfile synPos
    butProfile = zeros(length(isAx),maxRad);
    
    tic
    for s = 1:length(isAx);
        targSyn = isAx(s);
        preDists = sqrt((subs(:,1)-synAnchors(targSyn,1)).^2 + (subs(:,2)-synAnchors(targSyn,2)).^2 + ...
            (subs(:,3)-synAnchors(targSyn,3)).^2);
        minDist = min(preDists);
        targPre = find(preDists == minDist,1);
        
        nearVox = find(preDists<=checkRad);
        if ~isempty(nearVox)
            %             sc = sc + 1;
            
            targPre = find(nearVox == targPre,1);
            nearDist = preDists(nearVox);
            nearSubs = subs(nearVox,:);
            minVox = min(nearSubs,[],1);
            nearSubs = round(nearSubs - repmat(minVox,[size(nearSubs,1) 1]))+1;
            maxVox = max(nearSubs,[],1);
            nearVol = zeros(maxVox);
            nearInd = sub2ind(maxVox,nearSubs(:,1),nearSubs(:,2),nearSubs(:,3));
            nearVol(nearInd) = 1;
            %image(squeeze(sum(nearVol,1))*10)
            
            %             synPos(sc,:) = synAnchors(targSyn,:);
            %             synCon(sc,:) = edges(targSyn);
            
            
            for b = 1:maxRad;
                filtVol = nearVol;
                ball = balls{b};
                filtVol(nearInd(nearDist>(b*2))) = 0;
                %Ic = convn(filtVol,ball,'same');
                Ic = fastCon(filtVol,ball);
                
                
                butProfile(s,b) = max(Ic(:));
                
                if (b>3) &(butProfile(s,b)<.1)
                    break
                end
                
                maxIc = max(Ic,[],3);
                
                if 1
                subplot(3,1,1)
                image(sum(nearVol,3)*10);
                subplot(3,1,2)
                image(sum(ball,3)*256/max(ball(:))/b/2)
%                 ylim([1 size(nearVol,1)])
%                 xlim([1 size(nearVol,2)])
                subplot(3,1,3)
                
                image(maxIc*256)
                pause(.01)
                end
            end
                       
        else
            'no near'
            butProfile(s,maxRad) = 0;
        end
                
    end % run all synapses
    toc
    axProfile{a} = butProfile;
    ax2syn{a} = isAx;
    %         mapSyn{a} = synPos;
    %         mapCon{a} = synCon;
    
end

butSize.axProfile = axProfile;

%% Analyize profiles
fitThresh = 0.5; %define fit to ball
butSize.fitThresh = fitThresh;

ballVol = zeros(length(balls),1);
for b = 1:length(balls)
    ball = balls{b};
    ballVol(b) = sum(ball(:)>0);
end

histButVolBin = [0:100:3000];
histButVol = zeros(length(axProfile),length(histButVolBin));
clear butVols
for a = 1:length(axProfile)
    a
    prof = axProfile{a};
    plot(prof'),pause(.01)
    goodProf = 0;
    clear butVol
    for p = 1:size(prof,1);
        checkProf = prof(p,:);
        ballFit = max(find(checkProf> fitThresh));
        
        if ~isempty(ballFit)
            butVol(p) = ballVol(ballFit)*checkProf(ballFit);
        else
            butVol(p) = 0;
        end
    end
    butVols{a} = butVol;
    histButVol(a,:) = hist(butVol,histButVolBin);
    
end

butSize.histButVolBin = histButVolBin;
butSize.butVols = butVols;
butSize.histButVol = histButVol;

subplot(2,1,1)
image(histButVol*4)
subplot(2,1,2)


%%
seedNum = length(seedList)
con = useList.con;
for s = 1:seedNum
    targ = find(useList.postList == seedList(s))
    getHist = histButVol(con(:,targ)>0,:);
    subplot(seedNum,1,s)
    bar(histButVolBin,sum(getHist,1))
end

%%
clear showRes
for a = 1:length(axList)
    
    showRes(a,1) = axList(a);
    showRes(a,2) = mean(butVols{a});
    showRes(a,3) = length(butVols{a});
end


for i = 1:length(butSize.butVols)
    
    butVols = butSize.butVols{i};
    butDiam = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2;
   
    butSize.butDiam{i} = butDiam;
end

save([MPN 'butSize4.mat'],'butSize')






