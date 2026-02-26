

seedCells = [108 201];


%% find cell prefs without exclustions
        for s = 1:length(seedCells)
             targ = find(cellList==seedCells(s));
             simMat = sqrt(con .* repmat(con(:,targ),[1 size(con,2)]));
             cellPref(s,:) = sum(simMat,1);
        end
        bar(cellPref')

%% Calculate what the preference of each post would be, absent a particular pre

postPrefMinusPre = zeros(size(con,1),size(con,2),2);
for i = 1: length(axList)
    usePre = setdiff(1:length(axList),i);
    sampCon = con(usePre,:);
    for s = 1:length(seedCells)
        targ = find(cellList==seedCells(s));
        simMat = sqrt(sampCon .* repmat(sampCon(:,targ),[1 size(con,2)]));
        postPrefMinusPre(i,:,s) = sum(simMat,1);
    end
end     


subplot(2,2,1)
image(postPrefMinusPre(:,:,1));
subplot(2,2,2)
image(postPrefMinusPre(:,:,2));
        
%% find ax prefrence given exclusion of post cell
clear sumMat
prePrefFromPost = zeros(size(con,1),size(con,2),2);
%normCellPref = cellPref./repmat(sum(cellPref,1),[2 1]);
for i = 1:length(axList)
    for p = 1:length(cellList)
        usePost = setdiff(1:length(cellList),p);
        sampPref = squeeze(postPrefMinusPre(i,usePost,:));
        sampCon = con(i,usePost);
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
isPrefInfo = ~isnan(postConPref) & ~isnan(preConPref);

minPref = min(max(postPrefMinusPre,[],3),max(prePrefFromPost,[],3));
isPrefInfo = isPrefInfo & (minPref>=minInfo);

subplot(2,2,3)
image(postConPref*100);
subplot(2,2,4)
image(preConPref*100);


%%
set(gcf,'Position',[100 100 950 900])

seedOne = 1:prod(size(con));
seedTwo = prod(size(con))+1:2*prod(size(con));
scatter(prePrefFromPost(seedOne),prePrefFromPost(seedTwo))
clf
hold on
for i = 1:size(con,1)
    for p = 1:size(con,2)
        if con(i,p) & isPrefInfo(i,p)
            scatter(postConPref(i,p,1),preConPref(i,p,1),con(i,p)*30,'k','filled');
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
for i = 1:size(con,1)
    for p = 1:size(con,2)
        if con(i,p) & isPrefInfo(i,p)
            %scatter(postConPref(i,p,1),preConPref(i,p,1),con(i,p)*30,'k');
        elseif isPrefInfo(i,p)
            scatter(postConPref(i,p,1),preConPref(i,p,1),10,'r','filled');
        end
    end
end

xlim([-.2 1.2])
ylim([-.2 1.2])


%%
conBin = [0:max(con(:))+1];
sameSeedOne = con((postConPref < .5) & (preConPref < .5) & isPrefInfo);
sameSeedTwo = con((postConPref > .5) & (preConPref > .5) & isPrefInfo);
difSeedOne = con((postConPref > .5) & (preConPref < .5) & isPrefInfo);
difSeedTwo = con((postConPref < .5) & (preConPref > .5) & isPrefInfo);
sameSeeds = [sameSeedOne; sameSeedTwo];
difSeeds = [difSeedOne; difSeedTwo];

histSame = histc(sameSeeds,conBin)/length(sameSeeds);
histDif = histc(difSeeds,conBin)/length(difSeeds);

bar(conBin,[histDif histSame])








