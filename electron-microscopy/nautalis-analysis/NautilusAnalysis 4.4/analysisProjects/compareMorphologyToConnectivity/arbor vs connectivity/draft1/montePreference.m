%%Analyze predictions of synapses based on cell preferences

clear all
clf
load('MPN.mat')
%[skelFile TPN] = uigetfile(MPN);
%load([TPN skelFile]);   %'skelOverlapPred'
load([MPN 'skelOverlapPred.mat']);   %'skelOverlapPred'




%%
binRes = 1; %resolution of bins used to predict overlap

useList = skelOverlapPred.useList;
cellList = useList.postList;
axList = useList.preList;
seedList = useList.seedList;
seedCells = [108 201];
con = skelOverlapPred.con;
sourceCon = con;
%%

for s = 1:length(seedCells)
    seedTarg(s) = find(cellList==seedCells(s));
    seedCon(:,s) = sourceCon(:,seedTarg(s));
end
% 
% useAx = sum(seedCon>0,2)>0;
% seedCon = seedCon(useAx,:);
% con = sourceCon(useAx,:);
% usePost = setdiff(1:size(con,2),seedTarg);
% con = con(:,usePost);

useCon = repmat(sum(seedCon>0,2),[1 size(sourceCon,2)]) & ...
    repmat(sum(sourceCon>0,1)>=2,[size(sourceCon,1) 1]);

binRes = 5;
predSyn = overlapPredictSynOfSubset(skelOverlapPred,binRes,useCon);


realPrefStat = checkPreference(sourceCon,predSyn,seedTarg,useCon,0)
testRand = checkPreference(sourceCon,predSyn,seedTarg,useCon,1)


realMean = realPrefStat.medPref;
realDif = realMean(1) - realMean(2);

if 0
%% randomize axons in front end
reps = 100
randMean = zeros(reps,2);
for r = 1:reps
    
    prefStat = checkPreference(sourceCon,predSyn,seedTarg,useCon,1);
    randMean(r,:) = prefStat.medPref;
    
end
randDif = randMean(:,1) - randMean(:,2);

histDifBin = [-1:.1:1];
histRandDif = hist(randDif,histDifBin);
bar([-1:.1:1],histRandDif/reps);
hold on
scatter(realDif,0,'filled','r')
hold off

rangeX(randDif)
median(randDif)
realDif
P = sum(randDif>=realDif)/reps

end

if 1
    %% Randomize by bootstrapping pair results
    reps = 10000;   
    sameSeedPref = realPrefStat.sameSeedPref;
    difSeedPref = realPrefStat.difSeedPref;
    realMean = [mean(sameSeedPref) mean(difSeedPref)];
    realDif = realMean(1) - realMean(2);
    realFrac =    [sum(sameSeedPref>0)/length(sameSeedPref) ...
            sum(difSeedPref>0)/length(difSeedPref)];
    
    numSame = length(sameSeedPref);
    numDif = length(difSeedPref);
    allRes = [sameSeedPref; difSeedPref];
 
    
    randMean = zeros(reps,2);
    for r = 1:reps
       pick = randperm(length(allRes));
       pickSame = allRes(pick(1:numSame));
       pickDif = allRes(pick(numSame+1:end));
        randMean(r,:) = [mean(pickSame) mean(pickDif)];
        randFrac(r,:)= [sum(pickSame>0)/length(pickSame) ...
            sum(pickDif>0)/length(pickDif)];
    end
    
    randDif = randMean(:,1) - randMean(:,2);

histDifBin = [-1:.1:1];
histRandDif = hist(randDif,histDifBin);
bar([-1:.1:1],histRandDif/reps);
hold on
scatter(realDif,0,'filled','r')
hold off

range = rangeX(randDif)
medianDifference = median(randDif)
realDif
P = sum(randDif>=realDif)/reps
  numSame  
  numDif
end

