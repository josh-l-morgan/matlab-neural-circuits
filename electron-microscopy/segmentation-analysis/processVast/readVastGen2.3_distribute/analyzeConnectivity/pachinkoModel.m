
clear all

SPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\skel\'

reps = 10000;



load([SPN 'terSubs_Ds8_Ds4_Look10.mat'])
synMat = terSubs.synMat;
touchMat = terSubs.touchMat;





%% Make overlap lookup table

sumOverlap = sum(touchMat(:));
[matY matX] = size(touchMat);
linMat = touchMat(:);
sumOverlap = sum(linMat);
lookupPair = zeros(sumOverlap,1);
prev = 0;
for i = 1:length(linMat);
    toEnd = prev + linMat(i);
    lookupPair(prev+1:toEnd) = i;
    prev = toEnd;
end
    
%% Make linear prediction
numSyn = sum(synMat(:));
predMat = touchMat * numSyn/sumOverlap;

realSqrC = squaredCluster(synMat);


%% Run model

for r = 1: reps

    
makeSyn = synMat * 0;
pickRand = ceil(rand(numSyn,1)*sumOverlap);
pickedTouch = lookupPair(pickRand);
histSyn = hist(pickedTouch,[1:1:length(linMat)]);
makeSyn(:) = histSyn;

randSqrC(r) = squaredCluster(makeSyn);

%image(makeSyn*10),pause(.01)

end

%% 
divHist = (0:.05:10);
histRandSqrC = hist(randSqrC,divHist)
bar(divHist,histRandSqrC/max(histRandSqrC))
hold on
scatter(realSqrC,0.01,'r');
hold off
xlim([0 10])
