
%%Check for unassigned bipolar cells
%%Determine if polarity of Ca responses are consistent with other cells. Is
%%it SNR dependent.
%%See if changing masks changes bimodality or polarity spread. 
%%Consolidate redundant ROIs

%{

1) Karl finishes cell 14, bip identification, ROI correlation
2) Try filtering ROIs by redundancy and edge proximity
3) Check distribution of physiological responses against other
distributions for RANGE and BIMODALITY


%}



if 0
    
    datFold = [glob.datDir 'Analysis\Data\preproc\'];
    SPN =datFold;
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
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
testLengths = 0:1:100;
L = length(testLengths);

caPol = ROI.Polarity;
SNR = ROI.SNR;
SNR(isnan(SNR)) = 0;
SNRcol = colorProp(SNR,'STD');

roiCids = ptDat(:,3);
numRoi = size(ptDat,1);

allPred = zeros(size(ptDat,1),L,length(onScale),length(offScale),length(noise));
goodRoi = zeros(numRoi,1);
goodCid = zeros(length(SOI.cids),1);
useCid = SOI.cids;

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
            
            predPol = zeros(length(isCid),L,length(onScale),length(offScale));
            for c = 1:L
                lc = testLengths(c);
                
                W = exp(-d/lc); % Apply length constant
                for n = 1:length(noise)
                    for s = 1:length(onScale)
                        for s2 = 1:length(offScale)
                            Won = W .* onMat * onScale(s) + noise(n);
                            Woff = W .* offMat * offScale(s) + noise(n);
                            
                            sumOn = sum(Won,1) ;
                            sumOff = sum(Woff,1) ;
                            
                            predPol(:,c,s,s2,n) = (sumOff-sumOn)./(sumOff+sumOn);
                        end
                    end
                end
                
            end
            allPred(isCid,:,:,:,:) = predPol;
        end
    end
end




%% Find error
'finish reading out noise!!!!'
%useRoi = find(~isnan(sum(allPred,2)));

if filterBySNR
    goodRoi = goodRoi & (SNR>.07);
end
if filterByEdge
    goodRoi(roiNearEdge) = 0;
end
useRoi = find(goodRoi>0);


%%Set error weights
if weightErrors == 1
    ew = SNR;% * 0 + 1;
else %dont weight errors
    ew = SNR * 0 + 1;
end
ew = ew(goodRoi>0);
ew = ew/mean(ew);
ewM = repmat(ew,[1 size(usePred,2) size(usePred,3) size(usePred,4) size(usePred,5)]);


usePred = allPred(useRoi,:,:,:,:);
caPolN = caPol(useRoi);% - mean(caPol(useRoi));
if standardize
    caPolN = caPolN/std(caPolN);
end
usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
if standardize
    stdUsePred = std(usePredN,1);
    usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1 1]);
end
caPolMatN = repmat(caPolN,[1 L size(usePred,3) size(usePred,4) size(usePred,5)]);
difMatN = caPolMatN - usePredN;
meanErrN = mean(abs(difMatN));
rmseN = sqrt(mean(difMatN.^2,1));
rmseN = rmseN/ mean(abs(caPolN));

%%Get correlation coefficient
caPolMatM = caPolMatN - repmat(mean(caPolMatN,1),[size(caPolMatN,1) 1 1 1 1]);
usePredM = usePredN - repmat(mean(usePredN,1),[size(usePredN,1) 1 1 1 1]);
covN = mean(caPolMatM .* usePredM .* ewM,1);
std1 = std(caPolMatN,1);
std2 = std(usePredN,1);
cc = covN./(std1.*std2);

if matchType == 1 
    matchVal = (max(abs(rmseN(:))) - rmseN)./((std1 ));
else
    matchVal = cc;
end

%%Apply weights


minE = max(matchVal(:))
indE = find(matchVal==minE);
[by bx bs1 bs2 bn] = ind2sub(size(matchVal),indE);
bym = mean(by);
bxm = mean(bx);
bs1m = mean(bs1);
bs2m = mean(bs2);
bnm = mean(bn);


bestOnScale = onScale(round(bs1m))
bestOffScale = offScale(round(bs2m))
bestLength = testLengths(round(bxm))
bestNoise = noise(round(bnm))


clf
subplot(2,1,1)
showErr = matchVal(:,:,round(bs1m),round(bs2m),round(bnm));
plot(testLengths,showErr)
title(sprintf('best length = %0.2f, OnScale = %0.2f,OffScale = %0.2f, noise = %0.4f',...
    bestLength,bestOnScale,bestOffScale,bestNoise))

showPred = usePredN(:,:,round(bs1m),round(bs2m),round(bnm));
%showPred = squeeze(usePredN(:,round(bxm),:,round(bnm))); %% show all on scaling
%showPred = squeeze(usePredN(:,round(bxm),round(bzm),:)); %% show all noise



%%Show shift in matching with length constants
if 1
    subplot(2,2,3)
    scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    scat1.CData = SNRcol(useRoi,:);
    title(sprintf('length constant %0.2f',testLengths(2)))
    xlim([-1 1])
    ylim([-1 1])
    set(gca,'color','k')
    
    drawnow
    
    subplot(2,2,4)
    
    
    if 1
        for r = 2
            for s = 1:size(showPred,2)
                scat1 = scatter(showPred(:,s),caPolN,15,'filled');
                scat1.CData = SNRcol(useRoi,:);
                title(sprintf('length constant %0.2f',testLengths(s)))
                %                 xlim([-1 1])
                %                 ylim([-1 1])
                set(gca,'color','k')
                drawnow
                pause(.1)
            end
            for s = size(showPred,2):-1:1
                scat1 = scatter(showPred(:,s),caPolN,15,'filled');
                scat1.CData = SNRcol(useRoi,:);
                title(sprintf('length constant %0.2f',testLengths(s)))
                %                 xlim([-1 1])
                %                 ylim([-1 1])
                set(gca,'color','k')
                drawnow
                %pause(.1)
            end
            
        end
    end
    
end
return


%% Compare errror to SNR


bestErrors = difMatN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));

scatter(SNR(useRoi),bestErrors)

caN = caPolMatN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));
uPN = usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));

scatter3(caN,uPN,SNR(useRoi))



%% Check each roi
useRoi = find(goodRoi>0);
clf
for r = 1:length(useRoi);
    
    usePredR = usePredN(r,:,:,:,:);
    caPolNR = caPolN(r);% - mean(caPol(useRoi));
    dif = abs(usePredR-caPolNR);
        
    minE = min(dif(:));
    indE = find(dif==minE);
    [by bx bs1 bs2 bn] = ind2sub(size(dif),indE);
    bymR = mean(by);
    bxmR = mean(bx);
    bs1mR = mean(bs1);
    bs2mR = mean(bs2);
    bnmR = mean(bn);
        
    bestOnScale = onScale(round(bs1m));
    bestOffScale = offScale(round(bs2m));
    bestLengths(r) = testLengths(round(bxmR));
    bestNoise = noise(round(bnm));
    
    subplot(2,1,1)
    hold off
    scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    scat1.CData = SNRcol(useRoi,:);
    title(sprintf('length constant %0.2f',testLengths(s)))
    %                 xlim([-1 1])
    %                 ylim([-1 1])
    set(gca,'color','k')
    hold on
    scat2 = scatter(usePredN(r,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN(r),150,'filled');
        
    subplot(2,1,2)
    hold on
    plot(dif)
    drawnow
    
end

useBest = (bestLengths>testLengths(1)) & (bestLengths<testLengths(end));
useBestLengths = bestLengths(useBest);

percentLengthsUsed = mean(useBest)*100

subplot(2,1,1)
hold off
hist(useBestLengths,testLengths)

medBestLengthsR = median(useBestLengths)
sortBL = sort(useBestLengths,'ascend');
bl95 = [sortBL(round(length(sortBL) * .025))      sortBL(round(length(sortBL) * .975))]

return

%% Show each cell

useCid = SOI.cids(goodCid>0);
pCol  = {'r' 'g' 'b' 'k' 'm'}
clf
for v = 1:length(useCid)
    vCid = useCid(v);
    useRoi = find((goodRoi>0) & (roiCids==vCid));
    
    usePred = allPred(useRoi,:,:,:,:);
    caPolMat = repmat(caPol(useRoi),1, L, size(usePred,3), size(usePred,4), size(usePred,5));
    difMat = caPolMat - usePred;
    rms = sqrt(mean(difMat.^2,1));
    meanErr = mean(abs(difMat),1);
    
    
    %%Set error weights
    if weightErrors == 1
        ew = SNR;% * 0 + 1;
    else %dont weight errors
        ew = SNR * 0 + 1;
    end
    ew = ew(useRoi>0);
    ew = ew/mean(ew);
    ewM = repmat(ew,[1 size(usePred,2) size(usePred,3) size(usePred,4) size(usePred,5)]);
    
        
    caPolN = caPol(useRoi);% - mean(caPol(useRoi));
    if standardize
        caPolN = caPolN/std(caPolN);
    end
    usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
    if standardize
        stdUsePred = std(usePredN,1);
        usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1 1]);
    end
    caPolMatN = repmat(caPolN,[1 L size(usePred,3) size(usePred,4) size(usePred,5)]);
    difMatN = caPolMatN - usePredN;
    meanErrN = mean(abs(difMatN));
    rmseN = sqrt(mean(difMatN.^2,1));
    rmseN = rmseN/ mean(abs(caPolN));
    
    %%Get correlation coefficient
    caPolMatM = caPolMatN - repmat(mean(caPolMatN,1),[size(caPolMatN,1) 1 1 1 1]);
    usePredM = usePredN - repmat(mean(usePredN,1),[size(usePredN,1) 1 1 1 1]);
    ewM = repmat(ew,[1 size(usePredM,2) size(usePredM,3) size(usePredM,4) size(usePredM,5)]);
    covN = mean(caPolMatM .* usePredM .* ewM,1);
    std1 = std(caPolMatN,1);
    std2 = std(usePredN,1);
    cc = covN./(std1.*std2);
    
    if matchType == 1
        matchVal = (max(abs(rmseN(:))) - rmseN)./((std1 ));
    else
        matchVal = cc;
    end
    
    minE = max(matchVal(:))
    indE = find(matchVal==minE);
    [by bx bs1 bs2 bn] = ind2sub(size(matchVal),indE);
    bym = mean(by);
    bxm = mean(bx);
    bs1m = mean(bs1);
    bs2m = mean(bs2);
    bnm = mean(bn);
    
    bestOnScale = onScale(round(bs1m))
    bestOffScale = offScale(round(bs2m))
    bestLength = testLengths(round(bxm))
    bestNoise = noise(round(bnm))
    
    showErr = matchVal(:,:,round(bs1m),round(bs2m),round(bnm));
    showMeanErr = meanErrN(:,:,round(bs1m),round(bs2m),round(bnm));
    
    
    showPred = usePredN(:,:,round(bs1m),round(bs2m),round(bnm));
    
    %subplot(length(useCid),1,v)
    subplot(2,1,1)
    hold on
    
    plot(testLengths,showMeanErr,pCol{v})
    title(sprintf('best length = %0.2f, OnScale = %0.2f, OffScale = %0.2f, noise = %0.2f',bestLength,bestOnScale,bestOffScale,bestNoise))
    hold on
    subplot(2,1,2)
    hold on
    plot(testLengths,showErr,pCol{v})
    title('root mean square error')
    hold on
    
    ylim([0 1])
    
    %title(sprintf('%d rois for cell %d',length(useRoi),vCid))
    drawnow
end








