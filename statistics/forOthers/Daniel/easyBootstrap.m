%function[res] = bootBatch(dat)


%% fake data
%c1 = batch, c2 = group, c3 = value
dat = [1 1 1;1 1 2; 1 1 1 ; 1 1 1; 1 2 0; 1 2 0 ; 1 2 0; 2 1 1;...
    2 1 2; 2 1 2; 2 2 0; 2 2 0; 2 2 0; 3 1 2; 3 1 3; 3 1 1; 3 2 0]


%% set up model
useRank = 1;
reps = 100000;
useBatch = [1];
sumGroup = 1; %group to sum ranks

n = size(dat,1);
batches = unique(dat(:,1));
groups = unique(dat(:,2));
clear batch batchSize batchVal batchGroup batchRank
for bu = 1:length(useBatch);
    b = useBatch(bu);
    batch{b} = dat(dat(:,1)==b,:);
    batchSize(b) = size(batch{b},1);
    batchVal{b} = batch{b}(:,3);
    batchGroup{b} = batch{b}(:,2);
    [sortBatch batchIDX] = sort(batchVal{b},'descend');
    batchRank{b} = batchIDX/batchSize(b);
    
end

%%Choose value to compare (raw or rank)
if useRank
    useVal = batchRank;
else
    useVal = batchVal;
end

%% ranksum
for bu = 1:length(useBatch);
    b = useBatch(bu);
    g1 = batchVal{b}(batchGroup{b}==groups(1));
    g2 = batchVal{b}(batchGroup{b}==groups(2));
    r0 = ranksum(g1,g2);
    r1 = ranksum(g1,g2,'tail','right');
    r2 = ranksum(g1,g2,'tail','left');
    disp(sprintf('batch %d: twoTail = %.05f, right = %.05f, left = %.05f',b,r0,r1,r2)) 
end


%% measure reals
Val = zeros(length(useBatch),length(groups));
for bu = 1:length(useBatch)
    b = useBatch(bu);
    newGroup = batchGroup{b};
    for g = 1:length(groups)
        Val(b,g) = sum(useVal{b}(newGroup==g));
    end
end


%% measure randos
rVal = zeros(length(useBatch),length(groups),reps);
for r = 1:reps
    for bu = 1:length(useBatch)
        b = useBatch(bu);
        
        pick =  randperm(batchSize(b));
        newGroup = batchGroup{b}(pick);
        for g = 1:length(groups)
            
            rVal(b,g,r) = sum(useVal{b}(newGroup==g));
            
        end
    end
end


%% test statistic

%%Real vals

sumVal = sum(Val,1);
%tRat = sumVal(1)/sum(sumVal,2);
tRat = sumVal(1,sumGroup);

%%Rand Vals
rSumVal = sum(rVal,1);
%rTRat = squeeze(rSumVal(1,1,:)./sum(rSumVal,2));
rTRat = squeeze(rSumVal(1,sumGroup,:));

%%Stats
sortRand = sort(rTRat);
realVal = tRat
meanRand = mean(sortRand)
CI95 = [sortRand(round(reps*.025)) sortRand(round(reps*.975))]
pLower = mean(sortRand<=realVal)
pHigher = mean(sortRand>=realVal)

%Show distribution
clf
maxR = max([sortRand; realVal]);
minR = min([sortRand; realVal]);
hRange = [minR: (maxR-minR)/100 : maxR];
hRand = hist(sortRand,hRange);
bar(hRange,hRand/reps)
hold on
scatter(realVal,0,'o','r','filled')
hold off












