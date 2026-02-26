
if 0
    clear
    datFold = [glob.datDir 'Analysis\Data\preproc\'];
    SPN =datFold;
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
end

minBip = 10;
testLengths = 0:.1:30;
L = length(testLengths);

caPol = ROI.Polarity;
SNR = ROI.SNR;
SNR(isnan(SNR)) = 0;
SNRcol = colorProp(SNR,'STD');

roiCids = ptDat(:,3);
numRoi = size(ptDat,1);

allPred = zeros(size(ptDat,1),L);
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
            
            predPol = zeros(length(isCid),L);
            for c = 1:L
                lc = testLengths(c);
                
                W = exp(-d/lc); % Apply length constant
                
                Won = W .* onMat;
                Woff = W .* offMat;
                
                sumOn = sum(Won,1);
                sumOff = sum(Woff,1);
                
                predPol(:,c) = (sumOff-sumOn)./(sumOff+sumOn);
                
                if c == L
                    
                    tooOn = sum(predPol(:,c)==-1);
                    if sum(tooOn)
                    end
                    
                    
                end
                
            end
            allPred(isCid,:) = predPol;
        end
    end
end



showPred = allPred;
showPred(isnan(showPred)) = 0;
colormap jet(100)
image((showPred+1)*50);


%%  Compare predictions to actual

wate = SNR;

%useRoi = find(~isnan(sum(allPred,2)));
useRoi = find(goodRoi>0);


usePred = allPred(useRoi,:);
caPolMat = repmat(caPol(useRoi),[1 L]);
difMat = caPolMat - usePred;
rms = sqrt(mean(difMat.^2,1));
meanErr = mean(abs(difMat),1);
plot(rms)


caPolN = caPol(useRoi);% - mean(caPol(useRoi));
caPolN = caPolN/mean(abs(caPolN));
usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
stdUsePred = mean(abs(usePredN),1);%std(usePredN,1);
usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1]);
caPolMatN = repmat(caPolN,[1 L]);
difMatN = caPolMatN - usePredN;
meanErrN = mean(abs(difMatN));
rmseN = sqrt(mean(difMatN.^2,1));
rmseN = rmseN/ mean(abs(caPolN));
showErr = rmseN;

clf
subplot(4,1,1)
plot(testLengths,meanErrN)
%ylim([0 max(meanErrN)])

%%Show shift in matching with length constants
if 1
    subplot(4,1,2)                
    scat1 = scatter(usePredN(:,2),caPolN,15,'filled');
    scat1.CData = SNRcol(useRoi,:);
    title(sprintf('length constant %0.2f',testLengths(2)))
%     xlim([-1 1])
%     ylim([-1 1])    
    set(gca,'color','k')

    drawnow
    
    subplot(4,1,3)
    scat1 = scatter(usePredN(:,end),caPolN,15,'filled');
    scat1.CData = SNRcol(useRoi,:);
        set(gca,'color','k')

    title(sprintf('length constant %0.2f',testLengths(end)))
%     xlim([-1 1])
%     ylim([-1 1])
    drawnow
    subplot(4,1,4)
    
    
    for r = 2
        for s = 1:L
            scat1 = scatter(usePredN(:,s),caPolN,15,'filled');
            scat1.CData = SNRcol(useRoi,:);
            title(sprintf('length constant %0.2f',testLengths(s)))
%             xlim([-1 1])
%             ylim([-1 1])
            set(gca,'color','k')
            drawnow
        end
        for s = L:-1:2
            scat1 = scatter(usePredN(:,s),caPolN,15,'filled');
            scat1.CData = SNRcol(useRoi,:);
            title(sprintf('length constant %0.2f',testLengths(s)))
%             xlim([-1 1])
%             ylim([-1 1])
            set(gca,'color','k')
            drawnow
        end
        
    end
end

return

%% Show each cell

%%Reconsider normalization and weighting


useCid = [2 3 4];%SOI.vgcCids(goodCid>0);
pCol  = {'r' 'g' 'b'}
clf
for v = 1:length(useCid)
    vCid = useCid(v);
    useRoi = find((goodRoi>0) & (roiCids==vCid));
    useWate  = wate(useRoi);
    usePred = allPred(useRoi,:);
    caPolN = caPol(useRoi);% - mean(caPol(useRoi));
    caPolN = caPolN/mean(abs(caPolN));
    usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
    stdUsePred = mean(abs(usePredN),1);
    usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1])*.2;
    caPolMatN = repmat(caPolN,[1 L]);
    difMatN = caPolMatN - usePredN;
    difMatN = difMatN .* repmat(useWate,[1 size(difMatN,2)]);
    
    meanErrN = mean(abs(difMatN));
    rmseN = sqrt(mean(difMatN.^2,1));
    
    showErr = rmseN;
    
    
    %subplot(length(useCid),1,v)
    subplot(2,1,1)
    hold on
    plot(testLengths,meanErrN,pCol{v})
    title('mean error')
    hold off
    subplot(2,1,2)
    hold on
    plot(testLengths,rmseN,pCol{v})
    title('root mean square error')
    hold off
    
    %title(sprintf('%d rois for cell %d',length(useRoi),vCid))
    
end








