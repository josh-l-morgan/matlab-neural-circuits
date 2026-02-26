function[datProb cumProb] = getDatProbs(dat);



totalAx = sum(dat,2)
apOc = 1:length(totalAx);
totalAp = totalAx.*apOc';
sumAp = sum(totalAp)


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

sumSyn = sum(realSyn);



%% Find probabilityies
reps = 100000
fillFrac = sumSyn/sumAp;
synBin = [0:1:size(dat,2)*3];
clear datProb
for i = 1:length(apOc)
   modSyn = rand(reps,i)<=fillFrac;
   sumMod = sum(modSyn,2);
   modHist = hist(sumMod,synBin);
    datProb(i,:) = modHist/reps;
    
    
    
end
    
datProb
%% find cumulative probabilities
clear cumProb
for i = 1:size(datProb,2)
    
   cumProb(:,i) = sum(datProb(:,i:end),2);
end
cumProb




