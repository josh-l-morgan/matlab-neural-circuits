clear all

MPN = GetMyDir

SPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\matOut_AutoCB\';
load([SPN 'cbObj.mat'])

targCell = 108
maxLook = 50;

cbScale = [16 * 8 * .0046 16 * 8 * .004    16 * .030]
cbScale .* cbObj.ImageSize

%% load data
TPN = [MPN 'cbVSnetwork\'];
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])

anchorScale = [.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3) * 4];

for i = 1:length(dsObj)
    dsObj(i).subs = scaleSubs(double(dsObj(i).subs),voxelScale);
end

allAnchors = scaleSubs(obI.cell.anchors,anchorScale);
cellRefList = obI.cell.name;

%% Define variables
allSubs = cat(1,dsObj.subs);
maxSub = max(allSubs,[],1)+10;
clear allSubs

useList = obI2cellList_seedInput(obI,targCell);
preList = useList.preList;
postList = useList.postList;

cellIDs = obI.nameProps.cellNum(obI.cell.mainObID);
tcrList = cellIDs(obI.nameProps.tcr(obI.cell.mainObID));
postList = tcrList;



synMat = useList.con;

postAnchors = zeros(length(postList),3);
for i = 1:length(postList)
    postAnchors(i,:) = allAnchors(find(cellRefList == postList(i)),:);
end


%% find closest point and total length within 2x dist of closest point

showPlots = 1
cb2net.targCell  = targCell;
cb2net.preList = preList;
cb2net.postList = postList;



%%
        cellList = postList;
        cb2net.cellList = cellList;
        dsMaxSub = cbObj.ImageSize;
        for cellNum = 1: length(cellList)
            mid2 = getCellSubs(obI,dsObj,cellList(cellNum));
            cb2net.cells(cellNum).dsSubs = mid2;
        end
%}



%%
viewVol = zeros(cbObj.ImageSize,'uint8');
maxCBid = max(cbObj.vol(:));
for post = 1: length(postList)
    
    dsSubs = cb2net.cells(find(cb2net.cellList == postList(post))).dsSubs;
    dsSubs = fix(scaleSubs(dsSubs,1./cbScale))+1;
    dsSubs =  wallSubs(dsSubs,[1 1 1 ; cbObj.ImageSize]);
    dsInd = sub2ind(cbObj.ImageSize,dsSubs(:,1),dsSubs(:,2),dsSubs(:,3));
    
    subplot(2,1,1)
    cbVals = cbObj.vol(dsInd);
    cbVals = cbVals(cbVals>0);
    uCB = unique(cbVals);
    showOverlap = 0;
    if length(uCB)== 0
        'missing CB'
        postCB(post) = 0;
        showOverlap = 1;
    elseif length(uCB) == 1
        postCB(post) = uCB;
                showOverlap = 1;

    else
        histCB = hist(cbVals,uCB);
        postCB(post) = uCB(find(histCB == max(histCB),1));
                showOverlap = 1;

    end
    
    
    if showOverlap
        %%temp view combo
        hist(cbVals)
        
        viewVol = viewVol * 0;
        viewVol(dsInd) = 1;
        rangeZ = 1;
        middleZ = mode(dsSubs(:,3));
        
        sampCell = double(viewVol(:,:,max(1,middleZ-rangeZ):min(dsMaxSub(3),middleZ + rangeZ)));
        sampCB =  double(cbObj.vol(:,:,max(1,middleZ-rangeZ):min(dsMaxSub(3),middleZ + rangeZ))>0);
        
        colSamp(:,:,1) = sum(sampCell,3)*30;
        colSamp(:,:,2) = sum(sampCB,3)*10;
        colSamp(:,:,3) = sum(sampCB,3) * 10;
        
        subplot(2,1,2)
        image(uint8(colSamp))
        postCB(post)
        pause(.1)
    end
    %}
    
    % synMat(pre,post) = sum(synPre==preList(pre) & synPost == postList(post));
    
end
%synMat(pre,post) = length(uniqueInd);  % how many unique synapse positions are there
%}



%% assign cell to cb
cbID = zeros(1,maxCBid);
cbID(postCB(postCB>0)) = postList(postCB>0);

% 
% cbSynMat = zeros(length(preList),length(cbID));
% for i = 1:length(cbID)
%     
%     if cbID(i)>0
%         targ = find(postList == cbID(i));
%         cbSynMat(:,i) = synMat(:,targ);
%     end
% end



%% Find distances of CB to seed cell
targCB = find(cbID == targCell);

cbProps = regionprops(cbObj.vol,'Centroid');
cbCenters = cat(1,cbProps.Centroid);
cbCenters = scaleSubs(cbCenters,cbScale);
cbDists = sqrt((cbCenters(:,1)-cbCenters(targCB,1)).^2 + ...
    (cbCenters(:,2)-cbCenters(targCB,2)).^2 + ...
    (cbCenters(:,3)-cbCenters(targCB,3)).^2);


cbDat.cbID = cbID;
cbDat.postCB = postCB;
cbDat.cbCenters = cbCenters;
save([MPN 'cbDat.mat'],'cbDat');


%% find shared Syn
targInputs = cbSynMat(:,targCB);
cbSharedSyn = cbSynMat*0;
for i = 1:length(cbID);
    cbSharedSyn(:,i) = min(cbSynMat(:,i),targInputs);
end

cbSyn = sum(cbSharedSyn,1);
cbBridge = sum(cbSharedSyn>0,1);



%%

inNetwork = sum(cbSynMat,1)>0;

useCB = inNetwork &  (cbID ~= targCell);
subplot(2,2,1)
scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB),1),'.')
xlim([0 200])

subplot(2,2,3)
scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB)>0,1),'.')
xlim([0 200])


useCB = (cbID ~= targCell);
subplot(2,2,2)
scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB),1),'.')
xlim([0 200])

subplot(2,2,4)
scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB)>0,1),'.')
xlim([0 200])

%%
useCB = (cbID ~= targCell);
subplot(2,1,1)
scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB),1),'.')
xlim([0 200])
ylabel('Shared synapse number')
xlabel('Distance ove cell body to origin cell')

subplot(2,1,2)
scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB)>0,1),'.')
xlim([0 200])
xlabel('Distance ove cell body to origin cell')
ylabel('Shared axon number')



%% Display distance to synapse relationship

inNetwork = sum(cbSynMat,1)>0;

useCB = inNetwork &  (cbID ~= targCell);
subplot(2,2,1)
scatter(cbDists(useCB),sum(cbSynMat(:,useCB),1),'.')
xlim([0 200])

subplot(2,2,3)
scatter(cbDists(useCB),sum(cbSynMat(:,useCB)>0,1),'.')
xlim([0 200])


useCB = (cbID ~= targCell);
subplot(2,2,2)
scatter(cbDists(useCB),sum(cbSynMat(:,useCB),1),'.')
xlim([0 200])

subplot(2,2,4)
scatter(cbDists(useCB),sum(cbSynMat(:,useCB)>0,1),'.')

xlim([0 200])
%% Threshold look area


closeCB = find((cbDists'<maxLook) & (cbID ~= targCell));

closeSyn = sum(cbSharedSyn(:,closeCB),1);
closeBridge = sum(cbSynMat(:,closeCB)>0,1);


subplot(2,1,1)
hist(closeSyn,[0:1:20])
subplot(2,1,2)
hist(closeBridge,[0:1:8])


%% Draw volume
subplot(1,1,1)
closeVol = zeros(cbObj.ImageSize);
conVol = closeVol;
targVol = closeVol;
for i = 1:length(closeCB)
    
    conVol(cbObj.PixelIdxList{closeCB(i)}) = cbSyn(closeCB(i));
    closeVol(cbObj.PixelIdxList{closeCB(i)}) = 1;
    
end

targVol(cbObj.PixelIdxList{targCB}) = 1;

viewDim = 1;
maxConVol = squeeze(max(conVol,[],viewDim));
maxCloseVol = squeeze(sum(closeVol,viewDim));
sumTargVol = squeeze(sum(targVol,viewDim));
maxCloseVol = (maxCloseVol>0) * 20 + maxCloseVol * 2;
maxCloseVol(maxCloseVol>150) = 150;

colComb = maxCloseVol +maxConVol*10 + sumTargVol*30;
colComb(:,:,2) = maxCloseVol +maxConVol*100;
colComb(:,:,3) = maxCloseVol;

image(uint8(colComb))

%{
fileName = sprintf('conSphere_look%d_Dim%d.png',maxLook,viewDim)
imwrite(uint8(colComb),[TPN fileName])
%}


%%  Reference sphere

sumAllCB = squeeze(sum(cbObj.vol>0,viewDim));
sumAllCB = (sumAllCB>0) * 40 + sumAllCB * 2;
sumAllCB(sumAllCB>150) = 150;

colComb =  maxCloseVol*2 +maxConVol*10 + sumTargVol*30;
colComb(:,:,2) = maxCloseVol * 2 +maxConVol*100;
colComb(:,:,3) = maxCloseVol  + sumAllCB;

image(uint8(colComb))

%{
fileName = sprintf('conSphere_Reference_Dim%d.png',maxLook,viewDim)
imwrite(uint8(colComb),[TPN fileName])
%}





targCell
%%
cb2net.synMat = synMat;
%         cb2net.touchMat = touchMat;
% cb2net.histMat = histMat;
% cb2net.minTouch = minTouch;





%%

useVals = touchMat>=0;
Xdat = touchMat(useVals);
Ydat = synMat(useVals);
[rho pc] = corr(Xdat,Ydat)

scatter(Xdat,Ydat,'.','r')
hold on
maxX = max(Xdat);
X = [0 maxX];
[fit1 gof] = fit(Xdat,Ydat,'poly1');
Y = X* fit1.p1 + fit1.p2;
line(X ,Y);
hold off


pause(2)

fitDat.Xdat = Xdat;
fitDat.Ydat = Ydat;
fitDat.fit1 = fit1;
fitDat.gof = gof;
fitDat.rho = rho;
fitDat.corrP = pc;

cb2net.fitDat = fitDat;


%% Show syn vs length

colMat = synMat*256/max(synMat(:));
colMat(:,:,2) = touchMat * 256/max(touchMat(:));
colMat(:,:,3) = touchMat * 0;

%image(uint8(colMat))

%%


toc
