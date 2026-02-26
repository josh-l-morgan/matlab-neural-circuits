%%Analyze predictions of synapses based on cell preferences
clear all
MPN = GetMyDir
load([MPN 'skelOverlapPred.mat']);   %'skelOverlapPred'
load([MPN 'obI.mat'])

skelOverlapPred.measures
predSyn = skelOverlapPred.results{3}.predSyn;

seedCells = [108 201];
seedList = seedCells;
crossoverAxons = [2032	2033	2034	2035]
noSkel = [2014 1026]

onlyUseTraced = 0;


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
synPref = (con - predSyn)./predSyn;
minPredictedSynapseNumber = .01;

%% seed mat
seedMat = con*0;
for i = 1:length(seedCells)
   seedMat(:,find(cellList ==seedCells(i))) = 1; 
end


%% Select cells
muchTraced = [106 107 108 109 111 112 117 120 123 129 133 134 148 156 159 162 163 169 ...
    170 201 203 205 206 207 210 212 213 215 216 218];
useTraced = muchTraced; %setdiff(muchTraced,seedCells);
[wasTraced ia ib] = intersect(cellList,useTraced);
tracedPost = con*0;
tracedPost(:,ia) = 1;


if ~onlyUseTraced
tracedPost = con*0+1;
end




%% find cell prefs without exclusions
for s = 1:length(seedCells)
    targ = find(cellList==seedCells(s));
    simMat = sqrt(sourceCon .* repmat(sourceCon(:,targ),[1 size(sourceCon,2)]));
    cellPref(s,:) = sum(simMat,1);
    image(simMat*10)
end
bar(cellPref')

%% Calculate what the preference of each post would be, absent a particular pre

pullCrossover = [];
for i = 1:length(crossoverAxons)
   pullCrossover(i) = find(axList == crossoverAxons(i));
end


postPrefMinusPre = zeros(size(sourceCon,1),size(sourceCon,2),2);
for i = 1: length(axList)
   %usePre = setdiff(1:length(axList),[pullCrossover i]);

    usePre = setdiff(1:length(axList),[pullCrossover i]);
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


for i = 1:length(seedCells)
   seedPos(i) = find(cellList == seedCells(i));
end
   


clear sumMat
prePrefFromPost = zeros(size(sourceCon,1),size(sourceCon,2),2);
%normCellPref = cellPref./repmat(sum(cellPref,1),[2 1]);
for i = 1:length(axList)
    for p = 1:length(cellList)
        usePost = setdiff(1:length(cellList),[p seedPos]);
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


%% Make Cross mat
crossMat = con*0;
crossMat(pullCrossover,seedPos(end)) = 1;
image(crossMat*1000)

%%
minInfo = 1;
postConPref = postPrefMinusPre(:,:,1)./sum(postPrefMinusPre,3);
preConPref = prePrefFromPost(:,:,1)./sum(prePrefFromPost,3);
minPref = min(max(postPrefMinusPre,[],3),max(prePrefFromPost,[],3));

isPrefInfo = ~isnan(postConPref) & ~isnan(preConPref);
isPrefInfo = isPrefInfo & ~isnan(synPref) & ~isinf(synPref);
isPrefInfo = isPrefInfo & (minPref>=minInfo) & (predSyn>=minPredictedSynapseNumber);

subplot(2,2,3)
image(postConPref*100);
subplot(2,2,4)
image(preConPref*100);


%%
preConPref(crossMat>0)
con(crossMat>0)



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

hold on
for c = 1:length(pullCrossover)
    i = pullCrossover(c);
    for p = 1:size(sourceCon,2)
        if sourceCon(i,p) & isPrefInfo(i,p) & ~sum(seedPos == p)
            scatter(postConPref(i,p,1),preConPref(i,p,1),sourceCon(i,p)*30,'g','filled');
           
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
sameSeedOneMat = (postConPref < prefThresh) & (preConPref < prefThresh) & isPrefInfo & tracedPost;
sameSeedTwoMat = (postConPref > prefThresh) & (preConPref > prefThresh) & isPrefInfo & tracedPost;
difSeedOneMat = (postConPref > prefThresh) & (preConPref < prefThresh) & isPrefInfo & tracedPost;
difSeedTwoMat = (postConPref < prefThresh) & (preConPref > prefThresh) & isPrefInfo & tracedPost;

isSame = sameSeedOneMat | sameSeedTwoMat;
isDif = difSeedOneMat | difSeedTwoMat;
c
sameSeedCon = con(isSame);
difSeedCon = con(isDif);
sameSeedPred = predSyn(isSame);
difSeedPred = predSyn(isDif);
sameSeedPref = synPref(isSame);
difSeedPref = synPref(isDif);

%% bar preferences
clf



conBin = [min(synPref(isPrefInfo & tracedPost)): .5 :max(synPref(isPrefInfo & tracedPost))+1];

histSame = histc(sameSeedPref,conBin)/length(sameSeedPref) * 100;
histDif = histc(difSeedPref,conBin)/length(difSeedPref)*100;

bar(conBin(1:end),[histSame(1:end) histDif(1:end)],'barWidth',1)

mean(sameSeedPref),std(sameSeedPref)
mean(difSeedPref),std(difSeedPref)
ranksum(sameSeedPref,difSeedPref)

%% bar connectivity
clf


sameSeedOneMat = (postConPref < prefThresh) & (preConPref < prefThresh);
sameSeedTwoMat = (postConPref > prefThresh) & (preConPref > prefThresh) ;
difSeedOneMat = (postConPref > prefThresh) & (preConPref < prefThresh);
difSeedTwoMat = (postConPref < prefThresh) & (preConPref > prefThresh);

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




