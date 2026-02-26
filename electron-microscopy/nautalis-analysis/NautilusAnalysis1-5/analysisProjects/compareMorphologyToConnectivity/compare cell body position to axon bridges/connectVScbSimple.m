clear all

%MPN = GetMyDir
load('MPN.mat')

SPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\matOut_AutoCB\';
load([SPN 'cbObj.mat'])

targCell = 108

maxLook = 50;

%cbScale = [16 * 8 * .0046   16 * 8 * .004     16 * .030]

%% load data
TPN = [MPN 'cbVSnetwork\'];
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])


cbScale = obI.em.res .* [16 * 8    16 * 8      16] /1000
cbScale .* cbObj.ImageSize

%anchorScale = [.0184 0.016 0.030];
anchorScale = obI.em.res .* [ 4 4 1]/1000;
% voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3) * 4];
% 
% for i = 1:length(dsObj)
%     dsObj(i).subs = scaleSubs(double(dsObj(i).subs),voxelScale);
% end
% 

allAnchors = scaleSubs(obI.cell.anchors,anchorScale);
cellRefList = obI.cell.name;

%% Define variables
allSubs = cat(1,dsObj.subs);
maxSub = max(allSubs,[],1)+10;
clear allSubs

useList =  obI2cellList_seedInput_RGC_TCR(obI,targCell);
preList = useList.preList;
postList = useList.postList;
%postList = postList(postList<1000);
synMat = useList.con;
cellPref = seedPreferences(targCell,useList);

postAnchors = zeros(length(postList),3);
for i = 1:length(postList)
    postAnchors(i,:) = allAnchors(find(cellRefList == postList(i)),[1 2 3]);
end

targidx = find(postList == targCell);
targPos = postAnchors(targidx,:);

% %% Get median position
% clear postMedian postMode postMean
% for cellNum = 1: length(postList)
%     subs = getCellSubs(obI,dsObj,postList(cellNum));
%     subs = scaleSubs(subs,voxelScale);
% 
%     postMedian(cellNum,:) = median(subs,1);
%     postMode(cellNum,:) = mode(subs,1);
%     postMean(cellNum,:) = mean(subs,1);
% end


%% Get properties of linked cells
notSeed = setdiff(1:length(postList),targidx);
postDists = sqrt((postAnchors(notSeed,1)-targPos(1)).^2 + ...
    (postAnchors(notSeed,2)-targPos(2)).^2 + ...
    (postAnchors(notSeed,3)-targPos(3)).^2);
closePost = find(postDists <= maxLook);

binPost = [0:5:300];
histPostDists = histc(postDists,binPost)
postVol = 4/3 * pi * (binPost(2:end).^3);

plot(histPostDists(1:end-1)./postVol')


%% find shared Syn
postShareAx = cellPref.sharedAx(notSeed);
postShareSyn = cellPref.sharedSyn(notSeed);


%% Find distances of CB to seed cell

cbProps = regionprops(cbObj.vol,'Centroid');
cbCenters = cat(1,cbProps.Centroid);
cbCenters = scaleSubs(cbCenters(:,[2 1 3]),cbScale);

targCBpos =targPos

cbDists = sqrt((cbCenters(:,1)-targCBpos(1)).^2 + ...
    (cbCenters(:,2)-targCBpos(2)).^2 + ...
    (cbCenters(:,3)-targCBpos(3)).^2);

isClose = find(cbDists<=maxLook);
closeCBnum = length(isClose)-1;


%%
hieght = max(cbCenters(:,1));



scatter(cbCenters(:,2),hieght-cbCenters(:,1),'b','filled');
hold on
scatter(postAnchors(:,2),hieght-postAnchors(:,1),'c','filled')
scatter(cbCenters(isClose,2),hieght-cbCenters(isClose,1),'y','filled')
scatter(postAnchors(notSeed(closePost),2),hieght-postAnchors(notSeed(closePost),1),'r','filled')

scatter(targCBpos(:,2),hieght-targCBpos(:,1),'w','filled')
hold off

set(gca,'color','k')
%% Threshold look area

closeSyn = postShareSyn(closePost);
closeBridge = postShareAx(closePost);
nonSeedMat = synMat(:,notSeed);
closeMat = nonSeedMat(:,closePost);

if closeCBnum > length(closePost)
    closeSyn(closeCBnum) = 0;
    closeBridge(closeCBnum) = 0; 
    closeMat(1,closeCBnum) = 0; 
end


subplot(2,1,1)
hist(closeSyn,[0:1:20])
subplot(2,1,2)
hist(closeBridge,[0:1:8])

closeAx = sum(closeMat,2);
sum(closeAx>0)





'return'
return




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
closeMat = cbSynMat(:,closeCB);

subplot(2,1,1)
hist(closeSyn,[0:1:20])
subplot(2,1,2)
hist(closeBridge,[0:1:8])

closeAx = sum(closeMat,2);
sum(closeAx>0)


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

viewDim = 3;
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
