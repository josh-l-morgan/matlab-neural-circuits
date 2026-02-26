clear all

histBin = [0:1:5];

%x = synapse number 0-3
%y = apposition number 1-10
dat = [147 10 0 0 ;
    73 9 1 0;
    60 7 3 1;
    27 5 1 1;
    21 3 0 0;
    12 5 1 1;
    4 0 0 0;
    2 2 1 0;
    0 0 0 0;
    1 0 1 0];

 dat =     [
   150    13     0     0;
    74    15     0     0;
    66     7     1     0;
    27     4     2     0;
    24     4     2     2;
    11     5     0     0;
     7     2     1     1;
     4     0     0     0;
     2     2     1     0;
     2     0     1     0;
     4     0     0     0]
 
 
 [datProb cumProb]= getDatProbs(dat);



totalAx = sum(dat,2)
apOc = 1:length(totalAx);
totalAp = totalAx.*apOc'
sum(totalAp)


%% make appo list (ax id, synapse yes/no
ids = [];
realSyn = [];
ap = [];
synOc = [1 1 2 3];
synYes = [0 1 2 3];
axCount = 0;
for a = 1:size(dat,1);
    for s = 1:size(dat,2);
        ocNum = apOc(a);
        synProfile = zeros(1,ocNum);
        synProfile(1:synYes(s)) = 1;
        for i = 1:dat(a,s)
            startIn = length(ids);
            axCount = axCount+1;
            ids(startIn+1:startIn+ocNum) = axCount;
            realSyn(startIn+1:startIn+ocNum) = synProfile;
        end
    end
end
numAx = max(ids);

%% make Ap list
id2ap = hist(ids,1:max(ids));

%%  Count axons
realWithSyn = ids(realSyn>0);
realUAx = unique(realWithSyn);
clear realCountAx realCountAp
for i = 1:length(realUAx)
    realCountAx(i) = sum(realWithSyn==realUAx(i));
    realCountAp(i) = id2ap(realUAx(i));
end

realSynDist= hist(realWithSyn,[1:numAx]);
realApDist = hist(ids,[1:numAx]);

probIdx = sub2ind(size(cumProb),realApDist,realSynDist+1);
realCumProbs = cumProb(probIdx);
realMeanProb = mean(realCumProbs(realCumProbs<.5));

realMultSyn = sum(realCountAx-1);

realHistAx = hist(realCountAx,histBin)


%%  Model
reps = 10000;
modSyn = realSyn*0;
apNum = length(ids);
synNum = sum(realSyn);

histAx = zeros(reps,6);
allDat = zeros(size(dat,1),size(dat,2)*3,reps,'single');
for r = 1:reps
    newDat = zeros(size(dat,1),size(dat,2)*3);
    modSyn = modSyn*0;
    
    genRand = rand(1,apNum);
    [sortRand idx] = sort(genRand);
    pickRand = idx(1:synNum);
    modSyn(pickRand) = 1;
    trackSyn(r) = sum(modSyn);
    
    %modSyn = realSyn;
    
    withSyn = ids(modSyn>0);
   
synDist= hist(withSyn,[1:numAx]);
probIdx = sub2ind(size(newDat),realApDist,synDist+1);
useIdx = 1:max(probIdx);
histProbIdx = hist(probIdx,useIdx);
newDat(useIdx) = histProbIdx;

allDat(:,:,r) = newDat;

    
    if ~mod(r,1000)
        disp(sprintf('running %d of %d',r,reps))
    end
    
end

meanDat = mean(allDat,3);

%% Model second run to compare to meanDat

reps = 10000;
modSyn = realSyn*0;
apNum = length(ids);
synNum = sum(realSyn);

histAx = zeros(reps,6);
allDifDat = zeros(size(dat,1),size(dat,2)*3,reps,'single');
sumDifs = zeros(reps,1);
for r = 1:reps
    newDat = zeros(size(dat,1),size(dat,2)*3);
    modSyn = modSyn*0;
    
    genRand = rand(1,apNum);
    [sortRand idx] = sort(genRand);
    pickRand = idx(1:synNum);
    modSyn(pickRand) = 1;
    trackSyn(r) = sum(modSyn);
    
    %modSyn = realSyn;
    
    withSyn = ids(modSyn>0);
   
synDist= hist(withSyn,[1:numAx]);
probIdx = sub2ind(size(newDat),realApDist,synDist+1);
useIdx = 1:max(probIdx);
histProbIdx = hist(probIdx,useIdx);
newDat(useIdx) = histProbIdx;
difDat = (newDat-meanDat);
sumDifs(r) = sum(abs(difDat(:)))/2;

allDifDat(:,:,r) = difDat;

    
    if ~mod(r,1000)
        disp(sprintf('running %d of %d',r,reps))
    end
    
end



%%

bufDat = meanDat*0;
bufDat(1:size(dat,1),1:size(dat,2)) = dat;
realDifDat = bufDat-meanDat;
realSumDifs = sum(abs(realDifDat(:)))/2;

