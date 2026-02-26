%%Analyze predictions of synapses based on cell preferences

clear all
clf
MPN = GetMyDir
fileName = sprintf('%sshiftSkelOverlapPred.mat',MPN);
load(fileName); %'skelOverlapPred'

%%
binRes = 1; %resolution of bins used to predict overlap

useList = skelOverlapPred.useList;
cellList = useList.postList;
axList = useList.preList;
seedList = useList.seedList;
seedCells = seedList;
con = skelOverlapPred.con;
sourceCon = con;


%% Select cells
muchTraced = [106 107 108 109 111 112 117 120 123 129 133 134 148 156 159 162 163 169 ...
     170 201 203 205 206 207 210 212 213 215 216 218];
useTraced = muchTraced; %setdiff(muchTraced,seedCells);
[wasTraced ia ib] = intersect(cellList,useTraced);
tracedPost = con*0;
tracedPost(:,ia) = 1;
% 
% 
% difSynPred = overlapPredictSynOfSubset(skelOverlapPred,isDif);
% predSyn = skelOverlapPred.predSyn;
% synPref = (con - predSyn)./predSyn;
% minPredictedSynapseNumber = .01;


%% find cell prefs without exclusions
for s = 1:length(seedCells)
    targ = find(cellList==seedCells(s));
    simMat = sqrt(sourceCon .* repmat(sourceCon(:,targ),[1 size(sourceCon,2)]));
    cellPref(s,:) = sum(simMat,1);
    image(simMat*10)
end
bar(cellPref')

%% Calculate what the preference of each post would be, absent a particular pre

postPrefMinusPre = zeros(size(sourceCon,1),size(sourceCon,2),2);
for i = 1: length(axList)
    usePre = setdiff(1:length(axList),i);
    sampCon = sourceCon(usePre,:);
    for s = 1:length(seedCells)
        targ = find(cellList==seedCells(s));
        simMat = sqrt(sampCon .* repmat(sampCon(:,targ),[1 size(sourceCon,2)]));
        postPrefMinusPre(i,:,s) = sum(simMat,1);
    end
end

subplot(2,2,1)
image(postPrefMinusPre(:,:,1));
subplot(2,2,2)
image(postPrefMinusPre(:,:,2));

%% find ax prefrence given exclusion of post cell
clear sumMat
prePrefFromPost = zeros(size(sourceCon,1),size(sourceCon,2),2);
%normCellPref = cellPref./repmat(sum(cellPref,1),[2 1]);
for i = 1:length(axList)
    for p = 1:length(cellList)
        usePost = setdiff(1:length(cellList),p);
        sampPref = squeeze(postPrefMinusPre(i,usePost,:));
        sampCon = sourceCon(i,usePost);
        for s = 1:length(seedCells)
            simMat = sqrt(sampCon .* sampPref(:,s)');
            prePrefFromPost(i,p,s) = sum(simMat);
        end
    end
end

subplot(2,2,3)
image(prePrefFromPost(:,:,1));
subplot(2,2,4)
image(prePrefFromPost(:,:,2));


%% get preference filter
minInfo = 1;
postConPref = postPrefMinusPre(:,:,1)./sum(postPrefMinusPre,3);
preConPref = prePrefFromPost(:,:,1)./sum(prePrefFromPost,3);
minPref = min(max(postPrefMinusPre,[],3),max(prePrefFromPost,[],3));

isPrefInfo = ~isnan(postConPref) & ~isnan(preConPref);
%isPrefInfo = isPrefInfo & ~isnan(synPref) & ~isinf(synPref);
isPrefInfo = isPrefInfo & (minPref>=minInfo);% & (predSyn>=minPredictedSynapseNumber);

subplot(2,2,3)
image(postConPref*100);
subplot(2,2,4)
image(preConPref*100);


%% Draw bubbles
set(gcf,'Position',[100 100 950 900])

seedOne = 1:prod(size(sourceCon));
seedTwo = prod(size(sourceCon))+1:2*prod(size(sourceCon));
scatter(prePrefFromPost(seedOne),prePrefFromPost(seedTwo))
clf
hold on

for i = 1:size(sourceCon,1)
    for p = 1:size(sourceCon,2)
        if sourceCon(i,p) & isPrefInfo(i,p)
            scatter(postConPref(i,p,1),preConPref(i,p,1),sourceCon(i,p)*30,'k','filled');
        elseif isPrefInfo(i,p)
            % scatter(postConPref(i,p,1),preConPref(i,p,1),10,'r');
        end
        
        
    end
end

xlim([-.2 1.2])
ylim([-.2 1.2])
hold off
% %% Draw zeros for bubbles
% clf
% hold on
% for i = 1:size(sourceCon,1)
%     for p = 1:size(sourceCon,2)
%         if sourceCon(i,p) & isPrefInfo(i,p)
%             %scatter(postConPref(i,p,1),preConPref(i,p,1),sourceCon(i,p)*30,'k');
%         elseif isPrefInfo(i,p)
%             scatter(postConPref(i,p,1),preConPref(i,p,1),10,'r','filled');
%         end
%     end
% end
% 
% xlim([-.2 1.2])
% ylim([-.2 1.2])


%% Filter by favorite
%%Produce binary masks for con matrix that filters according to seed
%%preference

prefThresh = .5
sameSeedTwoMat = (postConPref < prefThresh) & (preConPref < prefThresh) & isPrefInfo;
sameSeedOneMat = (postConPref > prefThresh) & (preConPref > prefThresh) & isPrefInfo;
difSeedTwoMat = (postConPref > prefThresh) & (preConPref < prefThresh) & isPrefInfo;
difSeedOneMat = (postConPref < prefThresh) & (preConPref > prefThresh) & isPrefInfo;

isSame = sameSeedOneMat | sameSeedTwoMat;
isDif = difSeedOneMat | difSeedTwoMat;

sameSeedCon = con(isSame);
difSeedCon = con(isDif);

%% Collect masks

mask.isPrefInfo = isPrefInfo>0;
mask.tracePost = tracedPost>0;
mask.sameSeedOneMat = sameSeedOneMat>0;
mask.sameSeedTwoMat = sameSeedTwoMat>0;
mask.isSame = isSame>0;
mask.isDif = isDif>0;


%% bar connectivity
clf

isSame = sameSeedOneMat | sameSeedTwoMat;
isDif = difSeedOneMat | difSeedTwoMat;

sameSeedCon = con(isSame);
difSeedCon = con(isDif);

conBin = [0: 1 :max(con(:))+1];

histSame = histc(sameSeedCon,conBin)/length(sameSeedCon) * 100;
histDif = histc(difSeedCon,conBin)/length(difSeedCon)*100;

bar(conBin(2:end),[histSame(2:end) histDif(2:end)],'barWidth',1)
colormap jet
mean(sameSeedCon),std(sameSeedCon),std(sameSeedCon)/sqrt(length(sameSeedCon))
mean(difSeedCon),std(difSeedCon),std(difSeedCon)/sqrt(length(difSeedCon))
ranksum(sameSeedCon,difSeedCon)


%% Plot predictions
checkBins = [.1:.1:20];
for r = 1:length(checkBins);
    checkRes = checkBins(r);
    checkPred = overlapPredictSynOfSubset(skelOverlapPred,checkRes,(mask.isPrefInfo & mask.tracePost));
    checkSynError = conError(con,checkPred,(mask.isPrefInfo & mask.tracePost));
    recHitRate(r) = checkSynError.hitRate;
    recSynError(r) = checkSynError.synErrorPC;
        disp(sprintf('checking %d of %d, = %2.1f',r,length(checkBins),checkSynError.synErrorPC))

image(checkPred*20),pause(.1)
end

plot(checkBins,recHitRate,'r')
hold on
plot(checkBins,recSynError,'b')
hold off

%% Show cluster shifts

load([MPN 'obI.mat'])
conTo = makeConTo(obI,seedCells);
for s = 1:2;
    seed = conTo(s).targ;
    postList = setdiff(conTo(s).tcrList,seedCells);
    
    specFilt{s} = con* 0;
    for p = 1:length(postList)
        specFilt{s}(:,find(cellList == postList(p))) = 1;
    end
end
% 
% specFilt{1} = sameSeedOneMat;
% specFilt{2} = sameSeedTwoMat;
plotCol = {'r','b'};
clf

for f = 1:2
 %mask.isPrefInfo &
calcMask = mask.tracePost & specFilt{f} ;
synNumber = sum(con(calcMask));
shiftList = skelOverlapPred.shiftList;



useCon = con.*calcMask;
synapses = con2syn(useCon);
image(useCon*100),pause(1)

scrambleSyn(1:size(con,1),1:size(con,2),synapses);




clear recHitRate recSynError predCVRMSE
manySyn = 5;
for r = 1:length(shiftList);
    checkRes = 1;
    checkPred = overlapPredictSynOfSubset(skelOverlapPred,checkRes,(calcMask),'shift',r);
    numberPredicted(r) = sum(checkPred(calcMask));
%     predCVRMSE(r)  = cvrmse(checkPred(calcMask));
    checkPred = checkPred * synNumber/ numberPredicted(r);
    checkSynError = conError(con,checkPred,(calcMask));
    recHitRate(r) = checkSynError.hitRate;
    recSynError(r) = checkSynError.synErrorPC;
        disp(sprintf('checking %d of %d, = %2.1f',r,length(shiftList),checkSynError.synErrorPC))
    predCVRMSE(r) = checkSynError.predCVRMSE;
    synMany(r) = sum(checkPred(calcMask)>=manySyn);
%image(checkPred*20),pause(.1)
end

subplot(3,1,2)
plot(shiftList,predCVRMSE,plotCol{f})
ylim([0 5])
hold on
scatter(0,checkSynError.realCVRMSE,plotCol{f})

subplot(3,1,1)
plot(shiftList,numberPredicted,plotCol{f})
ylim([0 400])
hold on
scatter(0,checkSynError.realSynNumber,plotCol{f})



subplot(3,1,3)
plot(shiftList,synMany,plotCol{f})
ylim([0 40])
hold on
scatter(0,sum(con(calcMask)>=manySyn),plotCol{f})





end
hold off



%% Get predictions
binRes = 1;
predSyn = overlapPredictSynOfSubset(skelOverlapPred,binRes,(mask.isPrefInfo & mask.tracePost));
checkSynError = conError(con,predSyn,(mask.isPrefInfo & mask.tracePost))

image(synPred*20)

%%
sameSeedPred = predSyn(isSame);
difSeedPred = predSyn(isDif);
sameSeedPref = synPref(isSame);
difSeedPref = synPref(isDif);

%% Get prediction stats
predSyn = skelOverlapPred.predSyn;

predErr = conError(con,predSyn,(tracedPost & isPrefInfo))

sameSynPred = overlapPredictSynOfSubset(skelOverlapPred,isSame);
sameSynError1 = conError(con,predSyn,isSame)
sameSynError2 = conError(con,sameSynPred,isSame)

difSynPred = overlapPredictSynOfSubset(skelOverlapPred,isDif);

difSynError1 = conError(con,predSyn,isDif)
difSynError2 = conError(con,difSynPred,isDif)




%% bar preferences
clf


conBin = [min([sameSeedPref; difSeedPref]): .5 :max([sameSeedPref; difSeedPref])+1];
histSame = histc(sameSeedPref,conBin)/length(sameSeedPref) * 100;
histDif = histc(difSeedPref,conBin)/length(difSeedPref)*100;

bar(conBin(1:end),[histSame(1:end) histDif(1:end)],'barWidth',1)

mean(sameSeedPref),std(sameSeedPref)/sqrt(length(sameSeedPref))
mean(difSeedPref),std(difSeedPref)/sqrt(length(difSeedPref))
ranksum(sameSeedPref,difSeedPref)

return



%% Scatter pred to con

clf
conBin = [min(con(isPrefInfo & tracedPost)):max(con(isPrefInfo & tracedPost))+1];
plotMax = max(con(isPrefInfo))+1;

scatter(sameSeedPred,sameSeedCon,75,'b','filled')
hold on
scatter(difSeedPred,difSeedCon,75,'r','filled')
plot([0 plotMax],[0 plotMax],'k')

hold off

xlim([-.5 plotMax+.5])
ylim([-.5 plotMax+.5])

%%

cmap = jet(100);
cmap = cat(1,[0 0 0],cmap);
colormap(cmap)
image(con*10)
image(predSyn*10)




