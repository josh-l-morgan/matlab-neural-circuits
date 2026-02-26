
if 0
    
    datFold = [glob.datDir 'Analysis\Data\preproc\'];
    SPN =datFold;
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
end

standardize = 0;
onScale = [1];
offScale = [1];
noise = [0:.001:.2];

minBip = 10;
testLengths = 0:1:40;
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
                        Won = W .* onMat  + noise(n);
                        Woff = W .* offMat * onScale(s)+ noise(n);
                        
                        sumOn = sum(Won,1);
                        sumOff = sum(Woff,1);
                        
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
useRoi = find(goodRoi>0);

usePred = allPred(useRoi,:,:,:,:);
caPolMat = repmat(caPol(useRoi),1, L, size(usePred,3), size(usePred,4), size(usePred,5));
difMat = caPolMat - usePred;
rms = sqrt(mean(difMat.^2,1));
meanErr = mean(abs(difMat),1);


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

minE = min(rmseN(:))
indE = find(rmseN==minE);
[by bx bs1 bs2 bn] = ind2sub(size(rmseN),indE);
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
showErr = rmseN(:,:,round(bs1m),round(bs2m),round(bnm));
plot(testLengths,showErr)
title(sprintf('best length = %0.2f, OnScale = %0.2f, noise = %0.4f',bestLength,bestOnScale,bestNoise))

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
                pause(.1)
            end
            
        end
    
    
end

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
    
    minE = min(rmseN(:))
    indE = find(rmseN==minE);
    [by bx bs1 bs2 bn] = ind2sub(size(rmseN),indE);
    bym = mean(by);
    bxm = mean(bx);
    bs1m = mean(bs1);
    bs2m = mean(bs2);
    bnm = mean(bn);


bestOnScale = onScale(round(bs1m))
bestOffScale = offScale(round(bs2m))
bestLength = testLengths(round(bxm))
bestNoise = noise(round(bnm))


showErr = rmseN(:,:,round(bs1m),round(bs2m),round(bnm));
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

    ylim([0 3])
    
    %title(sprintf('%d rois for cell %d',length(useRoi),vCid))
    drawnow
end








