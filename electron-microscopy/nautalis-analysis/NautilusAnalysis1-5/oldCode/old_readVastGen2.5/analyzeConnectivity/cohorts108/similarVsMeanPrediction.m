
clear all

SPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\skel\'
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];

reps = 10000;



load([SPN 'terSubs_Ds8_Ds1_Look1.mat'])
synMat = terSubs.synMat;
touchMat = terSubs.touchMat;
touchMat = ones(size(touchMat));


zerMat = synMat;
zerMat(:,1) = 0;
axSum = sum(zerMat,2);
useAx = axSum>0;
zerMat = zerMat(useAx,:);
synMat = synMat(useAx,:);

axSum = sum(synMat,2);
simNormMat = zerMat./repmat(axSum,[1,size(zerMat,2)]);

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
recSim = zeros(size(synMat,1),size(synMat,1));
allPre = 1:size(synMat,1);
for pre = 1:size(normMat,1)
    otherPre = setdiff(allPre,pre);
    selfMat = repmat(simNormMat(pre,:),[size(normMat,1) 1]);
    %difMat = abs(normMat-selfMat);
    minMat = min(simNormMat,selfMat);
    minMat(:,1) = 0; %% Remove seed cell from similarity calculation
    simList = sum(abs(minMat),2);
    simOther = simList(otherPre).^3;
    normSimOther = simOther/sum(simOther);
    recSim(pre,otherPre) = normSimOther;
    for post = 1:size(normMat,2)
        predictions = normMat(otherPre,post);
        weightedPredictions = predictions.* normSimOther;
        simMat(pre,post) = sum(weightedPredictions);
    end
end

image(recSim*400)

% 
% 
% simMatPos = normMat * 0;
% for pre = 1:size(normMat,1)
%     otherPre = setdiff(allPre,pre);
%     selfMat = repmat(normMat(pre,:),[size(normMat,1) 1]);
%     difMat = abs(normMat-selfMat);
%     minMat = min(normMat,selfMat);
%     simList = sum(abs(minMat),2);
%     simOther = simList(otherPre);
%     for post = 1:size(normMat,2)
%         predictions = normMat(otherPre,post);
%         doesSynapse = predictions>0;
%         useSimOther = simOther(doesSynapse)
%         weightedPredictions = predictions.* normSimOther;
%         simMatPos(pre,post) = sum(weightedPredictions);
%     end
% end

%%


%%

%%


        Xdat = simMat(:,1)*100;
        Ydat = normMat(:,1)*100;
        [rho pc] = corr(Xdat,Ydat)
        
        scatter(Xdat,Ydat,'o','r','lineWidth',2)
        hold on
        maxX = max(Xdat);
        X = [0 maxX];
        [fit1 gof] = fit(Xdat,Ydat,'poly1');
        Y = X* fit1.p1 + fit1.p2;
        line(X ,Y);
        xlim([0 100])
        ylim([0 100])
        hold off
        
        
        pause(2)
        
        fitDat.Xdat = Xdat;
        fitDat.Ydat = Ydat;
        fitDat.fit1 = fit1;
        fitDat.gof = gof;
        fitDat.rho = rho;
        fitDat.corrP = pc;
        
        terSubs.fitDat = fitDat;
        




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
