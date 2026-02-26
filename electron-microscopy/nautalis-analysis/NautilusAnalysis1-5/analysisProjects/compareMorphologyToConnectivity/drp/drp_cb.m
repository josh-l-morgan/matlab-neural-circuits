
clear all
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'cbDat.mat']);
load([MPN 'obI.mat']);
seedList = [108 201];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);

%% Get distances to all cell bodies from seed cells
cbID = cbDat.cbID;
cbCenters = cbDat.cbCenters;
cbDists = zeros(length(cbID),length(seedList));
for s = 1:length(cbID)
    targCB = s;%find(cbID == seedList(s));
    
    cbDists(:,s) = sqrt((cbCenters(:,1)-cbCenters(targCB,1)).^2 + ...
        (cbCenters(:,2)-cbCenters(targCB,2)).^2 + ...
        (cbCenters(:,3)-cbCenters(targCB,3)).^2);
    
end

%% make volume Mask
maxCent = round(max(cbCenters,[],1));
minCent = max([1 1 1], round(min(cbCenters,[],1)));
bufMax = 10;
volMask = (ones(maxCent));
volMask(1:minCent(1)+bufMax,:,:) = 0;
volMask(maxCent(1)-bufMax:end,:,:) = 0;
volMask(:,1:minCent(2)+bufMax,:) = 0;
volMask(:,:,1:minCent(3)+bufMax) = 0;
volMask(:,maxCent(2)-bufMax:end,:) = 0;
volMask(:,:,maxCent(3)-bufMax:end) = 0;
maskInd = find(volMask>0);
[y x z] = ind2sub(size(volMask),maskInd);
maskSub = [y x z];


%% slow and flawed use of all cells
%{
binWidth = 10;
binRad = binWidth/2;
binRange = [binWidth/2:10:100];
binVol = zeros(length(cbID),length(binRange));
binCount = binVol;
binVox = binVol;
for s = 1:length(cbID)
    s
    cbPos = cbCenters(s,:);
    dists2Vox = sqrt((maskSub(:,1)-cbPos(1)).^2 + ...
        (maskSub(:,2)-cbPos(2)).^2 + ...
        (maskSub(:,3)-cbPos(3)).^2);
    
    
    dist2cb = sqrt((cbCenters(:,1)-cbPos(1)).^2 + ...
        (cbCenters(:,2)-cbPos(2)).^2 + ...
        (cbCenters(:,3)-cbPos(3)).^2);
    
    for i = 1:length(binRange)
        inSide = binRange(i)-binRad;
        outSide = binRange(i)+binRad;
        %         volumeSphere = ((4/3)*pi * (outSide/1000)^3) - ((4/3)*pi*(inSide/1000)^3);
        %         voxNumber = 0;
        binVox(s,i) =  sum((dists2Vox>inSide) + (dists2Vox<=outSide)==2);
        binCount(s,i) = sum((dist2cb>inSide) & (dist2cb<=outSide));
    end
end

binDensity = binCount(1:477,:)./binVox(1:477,:);
plot(mean(binDensity,1))
%}
%%

maxCheck = 50;
binWidth = 1;
binRad = binWidth/2;
binRange = [binWidth/2:1:maxCheck];
binVol = zeros(length(cbID),length(binRange));
binCount = binVol;
for s = 1:length(cbID)
        for i = 1:length(binRange)
        
        inSide = binRange(i)-binRad;
        outSide = binRange(i)+binRad;
        volumeSphere = ((4/3)*pi * (outSide/1000)^3) - ((4/3)*pi*(inSide/1000)^3);
        voxNumber = 0;
        currentDists = cbDists(:,s);
        binCount(s,i) = sum((currentDists>inSide) & (currentDists<=outSide));
        binVol(s,i) = volumeSphere;
    end
end



binDensity = binCount./binVol;
image(binDensity)

safeMin = repmat(min(cbCenters,[],1)+maxCheck,[length(cbID),1]);
safeMax = repmat(max(cbCenters,[],1)-maxCheck,[length(cbID),1]);

safeCB = sum((cbCenters>safeMin) & (cbCenters<safeMax),2) == 3;
safeMean = mean(binDensity(safeCB,:),1);
plot(binRange,safeMean)




%}