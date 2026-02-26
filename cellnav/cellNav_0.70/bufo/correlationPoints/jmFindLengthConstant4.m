 
if 0
clear 
SPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\individualPointMasksSpread\'

load([SPN 'ptDat.mat']);
load([SPN 'ROI.mat']);
load([SPN 'ROI2.mat']);
load([SPN 'SOI.mat']);
end

minBip = 10;
testLengths = 0:.1:30;
L = length(testLengths);

caPol = ROI2.Polarity;

roiCids = ptDat(:,3);
numRoi = size(ptDat,1);

allPred = zeros(size(ptDat,1),L);
goodRoi = zeros(numRoi,1);
goodCid = zeros(length(SOI.vgcCids),1);
useCid = SOI.vgcCids;

for v = 1:length(SOI.vgcCids)
    vCid = SOI.vgcCids(v);
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


%%

%useRoi = find(~isnan(sum(allPred,2)));
useRoi = find(goodRoi>0);

usePred = allPred(useRoi,:);
caPolMat = repmat(caPol(useRoi),[1 L]);
difMat = caPolMat - usePred;
rms = sqrt(mean(difMat.^2,1));
meanErr = mean(abs(difMat),1);
plot(rms)


caPolN = caPol(useRoi) - mean(caPol(useRoi));
caPolN = caPolN/std(caPolN)*.2;
usePredN = usePred - repmat(mean(usePred,1),[size(usePred,1) 1]);
stdUsePred = std(usePredN,1);
usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1])*.2;
caPolMatN = repmat(caPolN,[1 L]);
difMatN = caPolMatN - usePredN;
meanErrN = mean(abs(difMatN));
rmseN = sqrt(mean(difMatN.^2,1));
showErr = rmseN;

clf
subplot(2,1,1)
plot(testLengths,meanErrN)
subplot(2,1,2)


%%Show shift in matching with length constants
if 1
    for s = 1:L
        scatter(usePredN(:,s),caPolN);
        title(sprintf('length constant %0.2f',testLengths(s)))
        xlim([-1 1])
        ylim([-1 1])
        drawnow
    end
end


%% Show each cell

useCid = [2 3 4];%SOI.vgcCids(goodCid>0);
pCol  = {'r' 'g' 'b'}
clf
for v = 1:length(useCid)
    vCid = useCid(v);
    useRoi = find((goodRoi>0) & (roiCids==vCid));

usePred = allPred(useRoi,:);
caPolN = caPol(useRoi) - mean(caPol(useRoi));
caPolN = caPolN/std(caPolN)*.2;
usePredN = usePred - repmat(mean(usePred,1),[size(usePred,1) 1]);
stdUsePred = std(usePredN,1);
usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1])*.2;
caPolMatN = repmat(caPolN,[1 L]);
difMatN = caPolMatN - usePredN;
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







            
            