
%%Compare anatomically determined bipolar cell influence correlation 
%%Coefficient matches up with functional correlation coefficients for ROIs
%%Determine best length constant for matching
%%This code requires first running AutoCorrelation.mat


if 0
    
    datFold = [glob.datDir 'Analysis\Data\preproc\'];
    SPN =datFold;
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'COI.mat']);
    bipCids = COI.bipCids;
end

standardize = 0; %scale distributions by standard deviation before calculating error
matchType = 2; % 1 = rmse, 2 = cc
weightErrors = 1; % multiply covariance by weight (SNR) 
filterBySNR = 0;
filterByEdge = 1;


roiNearEdge = [66 67 68 69 71 120 144 70 198 ]; % manually identified by jm 1/6/2022
    
onScale = [1];
offScale = [1];
noise = [0];

minBip = 10;
testLengths = 1:1:100;
L = length(testLengths);

caPol = ROI.Polarity;
SNR = ROI.SNR;
SNR(isnan(SNR)) = 0;
SNRcol = colorProp(SNR,'STD');

roiCids = ptDat(:,3);
numRoi = size(ptDat,1);

allPred = zeros(length(bipCids),size(ptDat,1),L,length(noise));
goodRoi = zeros(numRoi,1);
goodCid = zeros(length(SOI.cids),1);
useCid = SOI.cids;


pos = SOI.pos;
allE = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2 + ...
    (pos(:,3)-pos(:,3)').^2);


%%Get weight
for v = 1:length(SOI.cids)
    vCid = SOI.cids(v);
    if sum(useCid==vCid)
        
        isCid = find(roiCids == vCid);
        useNodes = SOI.closeNode(isCid);
        if sum(useNodes)
            
            d = SOI.cell(v).d(:,useNodes);
            
            preSign = SOI.cell(v).preSign;
            numBip = sum(preSign>0);
            preSignMat = repmat(preSign,[1 size(d,2)]);
            onMat = preSignMat == 1;
            offMat = preSignMat == 2;
            if (sum(preSign==1)>=minBip) && (sum(preSign==2)>=minBip)
                fprintf('%d is good\n',vCid)
                goodRoi(isCid) = 1;
                goodCid(v) = 1;
            else
                fprintf('%d is bad\n',vCid)
            end
            
            pre = SOI.cell(v).syn.pre;
            clear bipSyn
            for b = 1:length(bipCids)
                bipSyn{b} = find(pre==bipCids(b));
            end

            bW = zeros(length(bipCids),size(d,2));
            predPol = zeros(size(bW,1),size(bW,2),L,length(noise));
            for c = 1:L
                lc = testLengths(c);
                for n = 1:length(noise)

                    W = exp(-d/lc)+ noise(n) ; % Apply length constant
                    for b = 1:length(bipCids)
                        bW(b,:) = sum(W(bipSyn{b},:),1);
                    end

                    predPol(:,:,c,n) = bW;

                end

            end
            allPred(:,isCid,:,:,:) = predPol;
        end
    end
end

bipR = zeros(numRoi,numRoi,L,length(noise));
for c = 1:L
    for n = 1:length(noise)
        bipR(:,:,c,n) = corrcoef(allPred(:,:,c,n));
    end
end

%% compare to autocorrelation



pos = SOI.pos;
allE = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2 + ...
    (pos(:,3)-pos(:,3)').^2);

allRoiInd = sub2ind([numRoi numRoi],allRoiY,allRoiX);
rMat = zeros(numRoi);
rMat(allRoiInd) = allR;

sameInd = allRoiInd(sameCell);
difInd = allRoiInd(difCell);

binW = .1;
hRange = [-.3:.01:1.1];
    difRoi = find(([1:numRoi] < [1:numRoi]')& (rMat~=0));


for c = 1:L
    for n = 1:length(noise)
        bR = bipR(:,:,c,n);
        bRSame = bR(sameInd);
        bRDif = bR(difInd);

        clear rBRSame  rBRDif
        for i = 1:length(hRange)
            l = hRange(i);
            rBRSame(i) = mean(rSame((bRSame>=(l-binW/2)) & (bRSame<(l+binW/2))));
            rBRDif(i) = mean(rDif((bRDif>=(l-binW/2)) & (bRDif<(l+binW/2))));
        end

        cc = corrcoef(bR(difRoi),rMat(difRoi));
        ccs(c,n) = cc(1,2);
                scatter(rMat(difRoi),bR(difRoi),5,'b','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)


        clf,
        subplot(2,1,1)
        hold on
        scatter(allE(difInd),bR(difInd),5,'b','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)
        scatter(allE(sameInd),bR(sameInd),5,'r','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)

        subplot(2,1,2),hold on

        %scatter(bR(difRoi),rMat(difRoi),5,'k','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)

        scatter(bRDif,rDif,5,'r','markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
        scatter(bRSame,rSame,5,'b','markeredgecolor',[0 0 1],'markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
        plot(hRange,rBRSame,'color',[0 0 0.4],'linewidth',2);
        plot(hRange,rBRDif,'color',[.4 0 0],'linewidth',2);
        title(sprintf('%d'),testLengths(c))
        drawnow
    end
end

clf
subplot(3,1,1)
plot(testLengths,ccs(:,1),'k')
maxCC = max(ccs(:));
[mt mn] = find(ccs==maxCC,1);

bestNoise = noise(mn);
bestLambda = testLengths(mt);
title(sprintf('the best match is with lambda = %0.1f, noise = %0.2f',bestLambda,bestNoise))


bR = bipR(:,:,mt,mn);
bRSame = bR(sameInd);
bRDif = bR(difInd);

clear rBRSame  rBRDif
for i = 1:length(hRange)
    l = hRange(i);
    rBRSame(i) = mean(rSame((bRSame>=(l-binW/2)) & (bRSame<(l+binW/2))));
    rBRDif(i) = mean(rDif((bRDif>=(l-binW/2)) & (bRDif<(l+binW/2))));
end

cc = corrcoef(bR(difRoi),rMat(difRoi));
ccs(c,n) = cc(1,2);

subplot(3,1,2)
hold on
scatter(allE(difInd),bR(difInd),5,'b','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)
scatter(allE(sameInd),bR(sameInd),5,'r','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)

subplot(3,1,3),hold on

%scatter(bR(difRoi),rMat(difRoi),5,'k','markerfacecolor','flat','markeredgealpha',.2,'MarkerFaceAlpha',0.1)

scatter(bRDif,rDif,5,'r','markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
scatter(bRSame,rSame,5,'b','markeredgecolor',[0 0 1],'markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
plot(hRange,rBRSame,'color',[0 0 0.4],'linewidth',2);
plot(hRange,rBRDif,'color',[.4 0 0],'linewidth',2);
title(sprintf('%d'),testLengths(c))
daspect([1 1 1])
drawnow


if 0
    %fDir = uigetdir;
    TPN = uigetdir
    filename = [TPN '\autoCor_EMtoCa'];
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])
    
end







