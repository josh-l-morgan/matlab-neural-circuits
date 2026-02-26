
clear all

SPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\skel\'
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];

reps = 10000;



load([SPN 'terSubs_Ds8_Ds1_Look1.mat'])
synMat = terSubs.synMat;
touchMat = terSubs.touchMat;
touchMat = ones(size(touchMat));





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
realSynDist = sort(synMat(:),'descend');

%% Run model

randSynDist = zeros(reps,length(synMat(:)));
for r = 1: reps

    
makeSyn = synMat * 0;
pickRand = ceil(rand(numSyn,1)*sumOverlap);
pickedTouch = lookupPair(pickRand);
histSyn = hist(pickedTouch,[1:1:length(linMat)]);
makeSyn(:) = histSyn;

randSqrC(r) = squaredCluster(makeSyn);
randSynDist(r,:) = sort(makeSyn(:),'descend');
%image(makeSyn*10),pause(.01)

end

meanRandSynDist = mean(randSynDist,1);
% for i = 1:size(randSynDist,1)
%    bar(randSynDist(i,:)),pause(.01) 
% end

bar([meanRandSynDist; realSynDist']')
plot([meanRandSynDist; realSynDist']','lineWidth',2)

%% 
divHist = (0:.05:10);
histRandSqrC = hist(randSqrC,divHist)
bar(divHist,histRandSqrC/max(histRandSqrC))
hold on
scatter(realSqrC,0.01,'r','LineWidth',5);
hold off
xlim([0 10])
