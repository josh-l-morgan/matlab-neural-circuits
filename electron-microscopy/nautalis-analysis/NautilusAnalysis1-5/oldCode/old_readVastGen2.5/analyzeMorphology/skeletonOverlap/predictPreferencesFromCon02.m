

seedCells = [108 201];
sourceCon = con;
finalCon = (con - predSyn)./predSyn;
minPredictedSynapseNumber = 1;

%% find cell prefs without exclustions
        for s = 1:length(seedCells)
             targ = find(cellList==seedCells(s));
             simMat = sqrt(sourceCon .* repmat(sourceCon(:,targ),[1 size(sourceCon,2)]));
             cellPref(s,:) = sum(simMat,1);
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



%%
minInfo = 1;
postConPref = postPrefMinusPre(:,:,1)./sum(postPrefMinusPre,3);
preConPref = prePrefFromPost(:,:,1)./sum(prePrefFromPost,3);
minPref = min(max(postPrefMinusPre,[],3),max(prePrefFromPost,[],3));

isPrefInfo = ~isnan(postConPref) & ~isnan(preConPref);
isPrefInfo = isPrefInfo & ~isnan(finalCon) & ~isinf(finalCon);
isPrefInfo = isPrefInfo & (minPref>=minInfo) & (predSyn>=minPredictedSynapseNumber);

subplot(2,2,3)
image(postConPref*100);
subplot(2,2,4)
image(preConPref*100);


%%
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
%%
clf
hold on
for i = 1:size(sourceCon,1)
    for p = 1:size(sourceCon,2)
        if sourceCon(i,p) & isPrefInfo(i,p)
            %scatter(postConPref(i,p,1),preConPref(i,p,1),sourceCon(i,p)*30,'k');
        elseif isPrefInfo(i,p)
            scatter(postConPref(i,p,1),preConPref(i,p,1),10,'r','filled');
        end
    end
end

xlim([-.2 1.2])
ylim([-.2 1.2])


%%
clf
conBin = [min(finalCon(isPrefInfo)): .5 :max(finalCon(isPrefInfo))+1];

sameSeedOne = finalCon((postConPref < .5) & (preConPref < .5) & isPrefInfo);
sameSeedTwo = finalCon((postConPref > .5) & (preConPref > .5) & isPrefInfo);
difSeedOne = finalCon((postConPref > .5) & (preConPref < .5) & isPrefInfo);
difSeedTwo = finalCon((postConPref < .5) & (preConPref > .5) & isPrefInfo);
sameSeeds = [sameSeedOne; sameSeedTwo];
difSeeds = [difSeedOne; difSeedTwo];

histSame = histc(sameSeeds,conBin)/length(sameSeeds);
histDif = histc(difSeeds,conBin)/length(difSeeds);

bar(conBin(1:end),[histDif(1:end) histSame(1:end)],'barWidth',1.5)

%%

clf
conBin = [min(finalCon(:)):max(finalCon(:))+1];

sameSeedOne2 = con((postConPref < .5) & (preConPref < .5) & isPrefInfo);
sameSeedTwo2 = con((postConPref > .5) & (preConPref > .5) & isPrefInfo);
difSeedOne2 = con((postConPref > .5) & (preConPref < .5) & isPrefInfo);
difSeedTwo2 = con((postConPref < .5) & (preConPref > .5) & isPrefInfo);
sameSeeds2 = [sameSeedOne2; sameSeedTwo2];
difSeeds2 = [difSeedOne2; difSeedTwo2];

scatter(sameSeeds,sameSeeds2,'r','filled')
hold on
scatter(difSeeds,difSeeds2,'g','filled')
plot([0 10],[0 10],'k')

hold off

xlim([0 10])
ylim([0 10])



