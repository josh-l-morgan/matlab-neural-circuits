%%Analyze predictions of synapses based on cell preferences
clear all
MPN = GetMyDir
load([MPN 'obI.mat'])

%% variables
seedList = [108 201]
crossoverAxons = [2032	2033	2034	2035]
noSkel = [2014 1026]

%% Get cells

useList = obI2cellList_seedInput(obI,seedList);
axList = useList.preList;
cellList = useList.postList;
synMat = useList.con;

skelOverlapPred.useList = useList;

synapses = obI.nameProps.edges;
edges = synapses(:,1:2);


axList = setdiff(axList,noSkel);
disp(sprintf('Results calculated without %d',noSkel));


%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end

sourceCon = con;
%% seed mat
seedMat = con*0;
for i = 1:length(seedList)
   seedMat(:,find(cellList ==seedList(i))) = 1; 
end

%% find cell prefs without exclusions
for s = 1:length(seedList)
    targ = find(cellList==seedList(s));
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
    for s = 1:length(seedList)
        targ = find(cellList==seedList(s));
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
        for s = 1:length(seedList)
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


%% Filter by favorite
%%Produce binary masks for con matrix that filters according to seed
%%preference

prefThresh = .5
sameSeedOneMat = (postConPref < prefThresh) & (preConPref < prefThresh) & isPrefInfo ;
sameSeedTwoMat = (postConPref > prefThresh) & (preConPref > prefThresh) & isPrefInfo ;
difSeedOneMat = (postConPref > prefThresh) & (preConPref < prefThresh) & isPrefInfo ;
difSeedTwoMat = (postConPref < prefThresh) & (preConPref > prefThresh) & isPrefInfo ;

isSame = sameSeedOneMat | sameSeedTwoMat;
isDif = difSeedOneMat | difSeedTwoMat;

sameSeedCon = con(isSame);
difSeedCon = con(isDif);


%% bar connectivity
clf


sameSeedOneMat = (postConPref < prefThresh) & (preConPref < prefThresh);
sameSeedTwoMat = (postConPref > prefThresh) & (preConPref > prefThresh) ;
difSeedOneMat = (postConPref > prefThresh) & (preConPref < prefThresh);
difSeedTwoMat = (postConPref < prefThresh) & (preConPref > prefThresh);

isSame = sameSeedOneMat | sameSeedTwoMat;
isDif = difSeedOneMat | difSeedTwoMat;


sameSeedCon = con(isSame & ~seedMat);
difSeedCon = con(isDif & ~seedMat);

conBin = [0: 1 :max(con(:))+1];

histSame = histc(sameSeedCon,conBin)/length(sameSeedCon) * 100;
histDif = histc(difSeedCon,conBin)/length(difSeedCon)*100;

bar(conBin(2:end),[histSame(2:end) histDif(2:end)],'barWidth',1)
colormap jet
mean(sameSeedCon),std(sameSeedCon),std(sameSeedCon)/sqrt(length(sameSeedCon))
mean(difSeedCon),std(difSeedCon),std(difSeedCon)/sqrt(length(difSeedCon))
ranksum(sameSeedCon,difSeedCon)


%% bar by axon

sameSeedConMat = con .* (isSame & ~seedMat);
difSeedConMat = con .*(isDif & ~seedMat);

sameSeedFrac = sum(sameSeedConMat,2)./ sum((isSame & ~seedMat),2);
difSeedFrac = sum(difSeedConMat,2)./ sum((isDif & ~seedMat),2);

mean(sameSeedFrac),std(sameSeedFrac)/sqrt(length(sameSeedFrac))
mean(difSeedFrac),std(difSeedFrac)/sqrt(length(difSeedFrac))
signrank(sameSeedFrac,difSeedFrac)

scatter(sameSeedFrac,difSeedFrac,5,1:length(difSeedFrac))
hold on
plot([0 3],[0 3])
hold off

xlim([0 3]);
ylim([0 3]);


