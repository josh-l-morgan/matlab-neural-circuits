
%%Compare anatomically determined bipolar cell influence correlation 
%%Coefficient matches up with functional correlation coefficients for ROIs
%%Determine best length constant for matching



%%
%%Find correlations between ROIs
%%What are correlations between nearby ROIs from different image planes

if 0
  figure
    global tis glob

    SPN = [glob.datDir 'Analysis\Data\preproc\'];
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'RawOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load([SPN 'COI.mat']);
    load([SPN 'GOI.mat']);PixResp
end


%% Get functional correlations
normalizeRunRs = 1;
useRank = 0;
gFilt = gausswin(1,5);
gFilt = gFilt/sum(gFilt);
pixExponent = 1; %Exponent to wich pixel values are raised during normalization
plot(gFilt)
useAllSameDif = 2;


%%Length constant variables
noise = 0;%[0:.03:.3];
minBip = 10;
testLengths = 1:1:50;
onlyUseHighCorr = 0; %if 1 select top onlyUseHighCorr of ROI pairs to use for finding length constant 
eucRange = [0 300]; % Range of distances between ROIs to allow use of



%sampRawBounds = [100 1500]; %%all
%sampRawBounds = [300 1500];
sampRawBounds = [100 300];  %% pre stimulus

usedRois = [];
allL = [];
allRnn = [];
allRnn = []; %not normalized
allE = [];
allry = [];
allrx = [];
allExp = [];
allR = [];
allRoiY = [];
allRoiX = [];
L = 0;

clf
for ry = 1:size(RawOI.useRois,1);
    for rx = 1:size(RawOI.useRois,2);

        pixResp = GOI.PixResp{ry,rx};
        useRois = GOI.useRois{ry,rx};
        rawPix = GOI.rawPixResp{ry,rx};
        roiCid = SOI.cid(useRois);
        uRoi = 1:length(useRois) < [1:length(useRois)]';
        [roiY roiX]  = find(uRoi);

        plot(rawPix')

        pix = rawPix(:,sampRawBounds(1):sampRawBounds(2))';
       subplot(1,2,1)
       image(pix*50/mean(pix(:)))
        pix = filter(gFilt,1,pix,[],1);
        pix = pix/mean(pix(:));
        pix = pix.^pixExponent;
        subplot(1,2,2)
       image(pix*50/mean(pix(:)))
           
       pause(.1)
       
        meanPix = mean(pix,1);
        pix = pix-repmat(meanPix,[size(pix,1) 1]);
        
        pixR = corrcoef(pix);

        r2rLin = SOI.r2rLin(useRois,useRois);
        r2rEuc = SOI.r2rEuc(useRois,useRois);
        
        N = sum(uRoi(:));
        uR = pixR(uRoi);
        allRnn(L+1:L+N) = uR;
        allL(L+1:L+N) = r2rLin(uRoi);
        allE(L+1:L+N) = r2rEuc(uRoi);
        allry(L+1:L+N) = ry;
        allrx(L+1:L+N) = rx;
        allRoiY(L+1:L+N) = useRois(roiY);
        allRoiX(L+1:L+N) = useRois(roiX);
        allExp(L+1:L+N) = ry*1000 + rx;
        L = L + N;
    end
end

%% normalize experiments
uExp = unique(allExp);

if normalizeRunRs
    %%Compare runs
    clf
    hold on
    expCol = hsv(length(uExp));

    oldR = allRnn;
    oldBounds = bounds95(allRnn,0.5)
    
    for i = 1:length(uExp)
        uExp(i)
        isExp = find(allExp==uExp(i));


        uR = oldR(isExp);
        uR = fitBounds(uR,oldBounds);
        allR(isExp) = uR;


        hold off
        scatter(allE,allRnn,4,'markerfacecolor',[0 0 0],...
            'markeredgecolor',[0 0 0],'markerfacecolor','flat',...
            'markeredgealpha',.5,'MarkerFaceAlpha',.2)

        hold on
        scatter(allE(isExp),allRnn(isExp),10,'markerfacecolor',[1 0 0],...
            'markeredgecolor',[1 0 0],'markerfacecolor','flat',...
            'markeredgealpha',1,'MarkerFaceAlpha',1)
        scatter(allE(isExp),allR(isExp),10,'markerfacecolor',[0 0 1],...
            'markeredgecolor',[0 0 1],'markerfacecolor','flat',...
            'markeredgealpha',1,'MarkerFaceAlpha',1)

        drawnow

    end
else
    allR = allRnn;
end
%% Display Functional correlations

sameCell = find(allL<inf);
difCell = find(allL==inf);

binW = 3;
hRange = [0:.1:150];


eSame = allE(sameCell);
rSame = allR(sameCell);
lSame = allL(sameCell);
eDif = allE(difCell);
lDif = allL(difCell);
rDif = allR(difCell);

clear reSame reDif rlDif rlSame reSameInt rlSameInt reDifInt
for i = 1:length(hRange)
    l = hRange(i);

    reSameVals = rSame((eSame>=(l-binW/2)) & (eSame<(l+binW/2)));
    reSame(i) = mean(reSameVals);
    reSameInt(i,:) = meanSem(reSameVals);
    rlSameVals = rSame((lSame>=(l-binW/2)) & (lSame<(l+binW/2)));
    rlSame(i) = mean(rlSameVals);
    rlSameInt(i,:) = meanSem(rlSameVals);
    reDifVals = rDif((eDif>=(l-binW/2)) & (eDif<(l+binW/2)));
    reDif(i) = mean(reDifVals);
    reDifInt(i,:) = meanSem(reDifVals);
    %rlDif(i) = mean(rDif((lDif>=(l-binW/2)) & (lDif<(l+binW/2))));
end

%%Make 3D surface relSame
binW2 = 5;
hRange2 =[0:1:150];

clear relSame
for i = 1:length(hRange2)
    l1 = hRange2(i);
    for o = 1:length(hRange2)
        l2 = hRange2(o);
        isE = (eSame>=(l1-binW2/2)) & (eSame<(l1+binW2/2));
        isL = (lSame>=(l2-binW2/2)) & (lSame<(l2+binW2/2));
        relSame(i,o) = mean(rSame(isE & isL));
        relSum(i,o) = sum(isE & isL);
    end
end
relSame(relSum<10) = nan;


%%Show data
clf
subplot(1,2,1)
hold on
scatter(eDif,rDif,2,'r','markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
scatter(lSame,rSame,2,'g','markeredgecolor',[0 .5 0],'markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
scatter(eSame,rSame,2,'b','markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
plot(hRange,reSame,'color',[0 0 0.4],'linewidth',2);
% plot(hRange,reSameInt(:,1),'color',[0 0 0.4],'linewidth',2);
% plot(hRange,reSameInt(:,3),'color',[0 0 0.2],'linewidth',2);
plot(hRange,reDif,'color',[.4 0 0],'linewidth',2);
plot(hRange,rlSame,'color',[0 0.4 0],'linewidth',2);

legend({' ',' ',' ','Euclidian distance between neurites of different cells',...
    'Euclidian distance between neurites of the same cell',...
    'Linear distance between neurites of the same cell'})
xlim([0 80])



subplot(1,2,2)
surf(hRange2,hRange2,relSame)
%daspect([1 1 .01])
daspect([1 1 100000])

view([90 -90])

pause(5)

if 0
    %fDir = uigetdir;
    TPN = uigetdir
    filename = [TPN '\autoCorr2_1-1500_2'];
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])
    
end


%% Find anatomical predictions
%%%%%%%%%%
standardize = 1; %scale distributions by standard deviation before calculating error
matchType = 2; % 1 = rmse, 2 = cc
weightErrors = 1; % multiply covariance by weight (SNR) 
filterBySNR = 0;
filterByEdge = 1;


roiNearEdge = [66 67 68 69 71 120 144 70 198 ]; % manually identified by jm 1/6/2022
    
L = length(testLengths);

caPol = GOI.Polarity;
SNR = GOI.SNR;
SNR(isnan(SNR)) = 0;
SNRcol = colorProp(SNR,'STD');

roiCids = GOI.roiCids;
numRoi = length(roiCids);

bipCids = COI.bipCids;
allPred = zeros(length(bipCids),numRoi,L,length(noise));
goodRoi = zeros(numRoi,1);
goodCid = zeros(length(GOI.cids),1);
useCid = GOI.cids;


pos = GOI.pos;
allE = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2 + ...
    (pos(:,3)-pos(:,3)').^2);


%%Get weight
for v = 1:length(GOI.cids)
    vCid = GOI.cids(v);
    if sum(useCid==vCid)
        
        isCid = find(roiCids == vCid);
        useNodes = GOI.closeNode(isCid);
        if sum(useNodes)
            
            d = SOI.cell(v).d(:,useNodes); %%GEt synapse to roi distances
            
%             preSign = SOI.cell(v).preSign;
%             numBip = sum(preSign>0);
%             preSignMat = repmat(preSign,[1 size(d,2)]);
%             onMat = preSignMat == 1;
%             offMat = preSignMat == 2;
%             if (sum(preSign==1)>=minBip) && (sum(preSign==2)>=minBip)
%                 fprintf('%d is good\n',vCid)
%                 goodRoi(isCid) = 1;
%                 goodCid(v) = 1;
%             else
%                 fprintf('%d is bad\n',vCid)
%             end
%             

            %%Generate list of synapses belonging to each possible bipolar
            %%cell in the dataset
            pre = SOI.cell(v).syn.pre;
            clear bipSyn
            for b = 1:length(bipCids)
                bipSyn{b} = find(pre==bipCids(b));
            end

            %Insert weight of synaptic influence of each bipolar cell into
            %each node
            bW = zeros(length(bipCids),size(d,2));
            predPol = zeros(size(bW,1),size(bW,2),L,length(noise));
            for c = 1:L
                lc = testLengths(c);
                for n = 1:length(noise)

                    W = exp(-d/lc)+ rand(size(d)) * noise(n) ; % make synapse to node weights using length constant
                    for b = 1:length(bipCids)
                        bW(b,:) = sum(W(bipSyn{b},:),1); %add together weights from each bipolar cell
                    end
                    predPol(:,:,c,n) = bW;

                end

            end
            allPred(:,isCid,:,:,:) = predPol;
        end
    end
end

meanPred = mean(allPred,5);
meanPred = mean(meanPred,4);
meanPred = mean(meanPred,3);
    



bipR = zeros(numRoi,numRoi,L,length(noise));

for c = 1:L
    for n = 1:length(noise)
        bW = allPred(:,:,c,n);
%         bWbounds = bounds95(bW(:));
%         bW(:) = fitBounds(bW(:));
        bipR(:,:,c,n) = corrcoef(bW);
    end
end

%% compare to autocorrelation



pos = SOI.pos;
allE = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2 + ...
    (pos(:,3)-pos(:,3)').^2);

allRoiInd = sub2ind([numRoi numRoi],allRoiY,allRoiX);
rMat = zeros(numRoi);
rMat(allRoiInd) = allR;
rMat(:) = fitBounds(rMat(:));
difRoi = find(([1:numRoi] < [1:numRoi]')& (rMat~=0));

if onlyUseHighCorr>0
    sMat = sort(rMat(rMat>0),'descend');
    matThresh = sMat(round(length(sMat)*onlyUseHighCorr));
    useRs = find(rMat>=matThresh);
end

useE = find((allE >= eucRange(1)) & (allE <= eucRange(2)));


sameInd = allRoiInd(sameCell);
difInd = allRoiInd(difCell);

binW = .1;
hRange = [-.3:.01:1];


clear ccs
for c = 1:L
    for n = 1:length(noise)
        bR = bipR(:,:,c,n);
        %bR(:) = fitBounds(bR(:));
        
        bRSame = bR(sameInd);
        bRDif = bR(difInd);

        %Define wich ROI pairs to use for correlation coeficient
        if useAllSameDif == 2
            roiForCC = sameInd;
        elseif useAllSameDif == 3;
            roiForCC = difInd;
        else
            roiForCC = [sameInd difInd];
        end


        if onlyUseHighCorr>0
            roiForCC = intersect(roiForCC, useRs);
        end

        roiForCC = intersect(roiForCC,useE);

        clear rBRSame  rBRDif
        for i = 1:length(hRange)
            l = hRange(i);
            rBRSame(i) = mean(rSame((bRSame>=(l-binW/2)) & (bRSame<(l+binW/2))));
            rBRDif(i) = mean(rDif((bRDif>=(l-binW/2)) & (bRDif<(l+binW/2))));
        end

        if ~useRank
            cc = corrcoef(bR(roiForCC),rMat(roiForCC));
            ccs(c,n) = cc(1,2);

        else
            urNum = length(roiForCC);
            [b bx] = sort(bR(roiForCC),'ascend');
            [r rx] = sort(rMat(roiForCC),'ascend');
            bRIdx(bx) = 1:urNum;
            rMatIdx(rx) = 1:urNum;
            cc = corrcoef(bRIdx,rMatIdx);
            avMax = urNum * urNum * 1/3; %?
            cc = (avMax-mean(abs(bRIdx-rMatIdx)))/avMax;
            ccs(c,n) = cc;
        end


        scatter(rMat(difRoi),bR(difRoi),5,'b','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)


        clf,
        subplot(2,1,1)
        hold on
        scatter(allE(difInd),bR(difInd),5,'b','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)
        scatter(allE(sameInd),bR(sameInd),5,'r','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)
       
        subplot(2,1,2),hold on

        %scatter(bR(difRoi),rMat(difRoi),5,'k','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)

        scatter(bRSame,rSame,5,'r','markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
        scatter(bRDif,rDif,5,'b','markeredgecolor',[0 0 1],'markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
        scatter(bR(roiForCC),rMat(roiForCC),15,'g','markeredgecolor',[0 1 0],'markerfacecolor','flat','markeredgealpha',1,'MarkerFaceAlpha',.2)

        plot(hRange,rBRSame,'color',[0 0 0.6],'linewidth',2);
        plot(hRange,rBRDif,'color',[.6 0 0],'linewidth',2);
        plot(hRange,hRange,'k')
        title(sprintf('%d'),testLengths(c))
        drawnow
    end
end

clf
subplot(2,2,1)
maxCC = max(ccs(:));
[mt mn] = find(ccs==maxCC,1);
plot(testLengths,ccs(:,mn),'k')

bestNoise = noise(mn);
bestLambda = testLengths(mt);
title(sprintf('the best match is with lambda = %0.1f, noise = %0.2f',bestLambda,bestNoise))


bR = bipR(:,:,mt,mn);
bRSame = bR(sameInd);
bRDif = bR(difInd);
bW = allPred(:,:,mt,mn);

clear rBRSame  rBRDif
for i = 1:length(hRange)
    l = hRange(i);
    rBRSame(i) = mean(rSame((bRSame>=(l-binW/2)) & (bRSame<(l+binW/2))));
    rBRDif(i) = mean(rDif((bRDif>=(l-binW/2)) & (bRDif<(l+binW/2))));
end


subplot(2,2,3)
hold on
scatter(allE(difInd),bR(difInd),5,'b','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)
scatter(allE(sameInd),bR(sameInd),5,'r','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)

subplot(1,2,2),hold on

%scatter(bR(difRoi),rMat(difRoi),5,'k','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)

scatter(bRDif,rDif,5,'b','markerfacecolor',[0 0 1],'markeredgealpha',.5,'MarkerFaceAlpha',.2)
scatter(bRSame,rSame,5,'r','markeredgecolor',[1 0 0],'markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
plot(hRange,rBRSame,'color',[0.4 0 0],'linewidth',2);
plot(hRange,rBRDif,'color',[0 0 0.4],'linewidth',2);
plot(hRange,hRange,'k')
title(sprintf('%d'),testLengths(c))
daspect([1 1 1])
drawnow
pause(1)

if 0
    %fDir = uigetdir;
    TPN = uigetdir
    filename = [TPN '\autoCor_EMtoCa'];
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])
    
end

return

%% Find interesting comparisons
getDifNum = 100;
%%Get most correlated ROIs on different cells
rDif = allR(difCell);
[sortDif idx] = sort(rDif,'descend');
topDifs = difCell(idx(1:100));

topPairs = [allRoiY(topDifs)' allRoiX(topDifs)'];


SPN = [glob.datDir 'Analysis\Data\preproc\'];
ROIMask = load([SPN 'maskDat.mat']);
ROIMask = ROIMask.maskDat;

%%What is the shared influence for these correlated ROIs?

topPtIlab = GOI.roiID(topPairs,1);
for i = 1:size(topPairs,1)
    
    clf

    subplot(3,1,1)
    bCC = bR(topPairs(i,1),topPairs(i,2))
    b1 = bW(:,topPairs(i,1));
    b2 = bW(:,topPairs(i,2));
    b1 = b1 * 256/max(b1);
    b2 = b2 * 256/max(b2);

    b12 = cat(2,b1,b2);
    image(b12')

    cmap = jet(256);
    cmap(1,:) = [0 0 0];
    colormap(gca,cmap)
    title(sprintf('cce = %0.2f',bCC))

    lab = GOI.roiID(topPairs(i,1),1);
    i1 = floor(lab/1000);
    i2 = mod(lab,1000);
    meanI1 = RawOI.meanIs{i1,i2} * .2;
    mask1 = (ROIMask(:,:,topPairs(i,1))>0)*30;
    colI1 = cat(3,meanI1,meanI1,meanI1);
    colI1(:,:,1) = colI1(:,:,1) + mask1;
    colI1(:,:,2) = colI1(:,:,2) - mask1;
    colI1(:,:,3) = colI1(:,:,3) - mask1;
    
    subplot(3,1,2)
    image(colI1)
    colormap(gca,gray(255))

    lab = GOI.roiID(topPairs(i,2),1);
    i1 = floor(lab/1000);
    i2 = mod(lab,1000);
    meanI2 = RawOI.meanIs{i1,i2} * .2;
    mask2 = (ROIMask(:,:,topPairs(i,2))>0)*30;
    colI2 = cat(3,meanI2,meanI2,meanI2);
    colI2(:,:,1) = colI2(:,:,1) + mask2;
    colI2(:,:,2) = colI2(:,:,2) - mask2;
    colI2(:,:,3) = colI2(:,:,3) - mask2;

    subplot(3,1,3)
    image(colI2)
    colormap(gca,gray(255))
    pause(1)

end





