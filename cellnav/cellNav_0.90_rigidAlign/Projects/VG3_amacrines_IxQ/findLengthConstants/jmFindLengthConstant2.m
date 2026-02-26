 


SPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\individualPointMasksSpread\'

load([SPN 'ptDat.mat']);
load([SPN 'ROI.mat']);
load([SPN 'ROI2.mat']);
load([SPN 'SOI.mat']);

testLengths = [1:50];
caPol = ROI2.Polarity;

roiCids = ptDat(:,3);

allPred = zeros(size(ptDat,1),length(testLengths));
for v = 1:length(SOI.vgcCids);
    isCid = find(roiCids == SOI.vgcCids(v));
    useNodes = SOI.closeNode(isCid);
    if sum(useNodes)
    
    d = SOI.cell(v).d(:,useNodes);
    preSign = SOI.cell(v).preSign;
    preSignMat = repmat(preSign,[1 size(d,2)]);
    onMat = preSignMat == 1;
    offMat = preSignMat == 2;
    
    predPol = zeros(length(isCid),length(testLengths));
    for c = 1:length(testLengths)
        lc = testLengths(c);
        
        W = exp(-d/lc); % Apply length constant
        
        Won = W .* onMat;
        Woff = W .* offMat;
        
        sumOn = sum(Won,1);
        sumOff = sum(Woff,1);
        
        predPol(:,c) = (sumOff-sumOn)./(sumOff+sumOn);
    end
    allPred(isCid,:) = predPol;
    end
    
end


showPred = allPred;
showPred(isnan(showPred)) = 0;
colormap jet(100)
image((showPred+1)*50);


%%

useRoi = find(~isnan(sum(allPred,2)));

usePred = allPred(useRoi,:);
caPolMat = repmat(caPol(useRoi),[1 length(testLengths)]);
difMat = caPolMat - usePred;
rms = sqrt(mean(difMat.^2,1));

[sortCal idx] = sort(caPolMat);


plot(rms)

for i = 1:numRoi
    bestF
    
    
end

            
            
            