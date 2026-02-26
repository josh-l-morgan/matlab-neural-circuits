
clear all

SPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\skel\'
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];

reps = 10000;



load([SPN 'terSubs_Ds8_Ds1_Look1.mat'])
synMat = terSubs.synMat;
touchMat = terSubs.touchMat;
touchMat = ones(size(touchMat));

axSum = sum(synMat,2);
normMat = synMat./repmat(axSum,[1,size(synMat,2)]);


%% 
meanMat = normMat * 0;
allPre = 1:size(normMat,1);
for pre = 1:size(normMat,1)
    otherPre = setdiff(allPre,pre);
    for post = 1:size(normMat,2)
        meanMat(pre,post) = mean(normMat(otherPre,post));
    end
end


simMat = normMat * 0;
for pre = 1:size(normMat,1)
    otherPre = setdiff(allPre,pre);
    selfMat = repmat(normMat(pre,:),[size(normMat,1) 1]);
    difMat = abs(normMat-selfMat);
    minMat = min(normMat,selfMat);
    simList = sum(abs(minMat),2);
    simOther = simList(otherPre);
    normSimOther = simOther/sum(simOther);
    for post = 1:size(normMat,2)
        predictions = normMat(otherPre,post);
        weightedPredictions = predictions.* normSimOther;
        simMat(pre,post) = sum(weightedPredictions);
    end
end

blankMat = normMat;
blankMat(:,1) = 0;
scatter(normMat(:),meanMat(:));

scatter(normMat(:),simMat(:));








% 
% 
% %% Make overlap lookup table
% 
% sumOverlap = sum(touchMat(:));
% [matY matX] = size(touchMat);
% linMat = touchMat(:);
% sumOverlap = sum(linMat);
% lookupPair = zeros(sumOverlap,1);
% prev = 0;
% for i = 1:length(linMat);
%     toEnd = prev + linMat(i);
%     lookupPair(prev+1:toEnd) = i;
%     prev = toEnd;
% end
%     
% %% Make linear prediction
% numSyn = sum(synMat(:));
% predMat = touchMat * numSyn/sumOverlap;
% 
% realSqrC = squaredCluster(synMat);
% realSynDist = sort(synMat(:),'descend');
% 
% %% Run model
% 
% randSynDist = zeros(reps,length(synMat(:)));
% for r = 1: reps
% 
%     
% makeSyn = synMat * 0;
% pickRand = ceil(rand(numSyn,1)*sumOverlap);
% pickedTouch = lookupPair(pickRand);
% histSyn = hist(pickedTouch,[1:1:length(linMat)]);
% makeSyn(:) = histSyn;
% 
% randSqrC(r) = squaredCluster(makeSyn);
% randSynDist(r,:) = sort(makeSyn(:),'descend');
% %image(makeSyn*10),pause(.01)
% 
% end
% 
% meanRandSynDist = mean(randSynDist,1);
% % for i = 1:size(randSynDist,1)
% %    bar(randSynDist(i,:)),pause(.01) 
% % end
% 
% bar([meanRandSynDist; realSynDist']')
% plot([meanRandSynDist; realSynDist']','lineWidth',2)
% 
% %% 
% divHist = (0:.05:10);
% histRandSqrC = hist(randSqrC,divHist)
% bar(divHist,histRandSqrC/max(histRandSqrC))
% hold on
% scatter(realSqrC,0.01,'r','LineWidth',5);
% hold off
% xlim([0 10])
