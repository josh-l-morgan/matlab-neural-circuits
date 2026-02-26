

clear all
MPN = 'D:\LGNs1\Analysis\segmentations\vastMip3\collectVastSubs_mat\';

mpSize = matlabpool('size');
if ~mpSize
    'opening matlab pool'
    matlabpool close force
    matlabpool
end

%%

% springRes = 'D:\LGNs1\Analysis\springDat\results\'
% saveName = 'randNoSeedsFilt1res_realNoSeedFilt1Wait30run25.mat'
% 
% load([springRes saveName])
results = allResults{1}

load('D:\LGNs1\Analysis\groupNodes\randSynRemoval\result_randSynRes_24min2.mat')

cladeNodes = cladeIt(results);
nodeIDs = results.nodeIDs;



%% generic
MPN ='D:\LGNs1\Analysis\segmentations\vastMip3\collectVastSubs_mat\';

%{
  [checkIDs checkProp]  = getList_cbArea(MPN);

        
    [checkIDs checkProp] = getList_mtaCount();
    
    [checkIDs checkProp] = getList_giantBoutonsCounts(MPN);
    [checkIDs checkProp] = getList_biconality;
    [checkIDs checkProp] = getList_giantBoutonsCounts(MPN);
[checkIDs checkProp]  = getList_axSeedCon();

%}

[checkIDs checkProp] = getList_mtaCount();



[checkIDs checkProp] = getList_giantBoutonsCounts(MPN);
useIDs = checkIDs<1000;
checkIDs = checkIDs(useIDs);
checkProp = checkProp(useIDs);
checkProp(checkProp>0) = 1;% checkProp(checkProp>0) + max(checkProp);


[checkIDs checkProp] = getList_pcLatRats;
checkProp = checkProp(:,1);

[checkIDs checkProp]  = getList_axSeedCon();


[overlap ai bi] = intersect(checkIDs, nodeIDs);
checkIDs = checkIDs(ai);
checkProp = checkProp(ai);

%% find clade differences


if (size(checkProp,1) == 1)
    checkProp = checkProp';
end


propRange = [min(checkProp(:)) max(sum(checkProp,2))];
propSpread =  max(sum(checkProp,2)) - min(checkProp(:));
histBin = [0:propSpread/30:propSpread];

%histBin = [0:.2:2];

meanProp = mean(checkProp,1);
difProp = checkProp - repmat(meanProp,[size(checkProp,1) 1]);
sumDif = sum(abs(difProp),2);
normDif = mean(sumDif);

%normProp = mean(sum(checkProp- repmat(mean(checkProp,1),[size(checkProp,1) 1]),2));
trackDif = trackNodeCladeDif(results,checkIDs,checkProp);
realTrackMean = trackDif.meanImprovement;


realDifRes = cladePropDif(cladeNodes,checkIDs,checkProp);
useVals = realDifRes.wcss;
useVals = useVals(realDifRes.allCounts>1);
realMeans = mean(useVals);
realHist = hist(useVals,histBin);

meanRealHist = mean(realHist,1);

clear randMeans randHist
reps = 100;
parfor i = 1:reps
    disp(sprintf('rand %d of %d',i,reps));
    randIDs = checkIDs(randperm(length(checkIDs)));
    %bar(checkProp(randperm(length(checkIDs)))*100)
    randDifRes =cladePropDif(cladeNodes,randIDs,checkProp);
    useVals = randDifRes.wcss;
    useVals = useVals(realDifRes.allCounts>1);
    randMeans(i) = mean(useVals);
    randHist(i,:) = hist(useVals,histBin);
    
    randTrackDif = trackNodeCladeDif(results,randIDs,checkProp);
    randTrackMean(i) = randTrackDif.meanImprovement;
    
    %     scatter(realDifRes.allCounts,realDifRes.sumDifs,'b')
    %     hold on
    %     scatter(randDifRes.allCounts, randDifRes.sumDifs,'r','.')
    %     hold off
    %     pause
    
    
end

meanRandHist = mean(randHist,1);

%%
subplot(3,1,1)
normHisBin = histBin / normDif;
bar(normHisBin,[meanRealHist/sum(meanRealHist); meanRandHist/sum(meanRandHist)]')


subplot(3,1,2)
histBin2 = [0:propSpread/100:propSpread];
histRandMeans = hist(randMeans,histBin2);
histRealMeans = hist(realMeans,histBin2);
normHistBin2 = histBin2 / normDif;
br = bar(normHistBin2,[ histRandMeans/sum(histRandMeans); histRealMeans/sum(histRealMeans)]')
set(br,'barwidth',1)

bar(normHistBin2,histRandMeans/sum(histRandMeans),'b')
hold on
bar(normHistBin2,[histRealMeans/sum(histRealMeans)],'r')
%set(br,'barwidth',1)
hold off



subplot(3,1,3)
biggestChange = max([realTrackMean randTrackMean]);
normHistBin3 = 0:biggestChange/10 : biggestChange;
histRealTrackMean = hist(realTrackMean,normHistBin3);
bar(normHistBin3, histRealTrackMean/max(histRealTrackMean),'r')
hold on
histRandTrackMean = hist(randTrackMean,normHistBin3);
bar(normHistBin3, histRandTrackMean/max(histRandTrackMean),'b')
hold off

%% stats
stat.P = sum(randMeans<=mean(realMeans))/length(randMeans);
stat.randMeans = randMeans;
stat.randHist = randHist;
stat.realMeans = realMeans;
stat.realHist = realHist;
stat.histBin = histBin;
stat








