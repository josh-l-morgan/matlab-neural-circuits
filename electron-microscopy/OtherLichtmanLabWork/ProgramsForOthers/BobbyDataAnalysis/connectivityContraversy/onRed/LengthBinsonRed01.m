
%% Set up variables
numRep = 1000;
maxSyn = 10;


%% Generate dummy data
axNum = 600;
axLength = rand(axNum,1)*3+1;
axOther = fix(rand(axNum,1)*6)-1;
axOther(axOther<0) = 0;
axRed = fix(rand(axNum,1)*6)-3;
axRed(axRed<0) = 0;

%% Analyze real data
realRat = axRed./(axRed+axOther);
realRat(isnan(realRat))=0;

realHistAxOther = hist(axOther,0:maxSyn);
realHistAxRed = hist(axRed,0:maxSyn);
realHistRat = hist(realRat,0:.1:1);

%%Count synapses
numOther = sum(axOther);
numRed = sum(axRed);
averageRat = numRed/(numRed+numOther);

%%Test Statistic
realDistFromAv = realRat-averageRat;
realStat =mean(realDistFromAv(realDistFromAv<0));

%%Display
subplot(4,2,1)
bar([0:1:maxSyn],realHistAxOther)
title('Real Other')
subplot(4,2,3)
bar([0:1:maxSyn],realHistAxRed)
title('Real Red')
subplot(4,2,5)
bar([0:.1:1],realHistRat)
title('Real Ratio')



%% Set up Test
sumLength = axLength*0;
for i = 1:axNum
    sumLength(i) = sum(axLength(1:i));
end
totalAx = sumLength(end);
sumLength = sumLength/totalAx;

%Create matrix where x is each synapse and y finds prop of axon

otherSumLength = repmat(sumLength,[1,numOther]);
redSumLength = repmat(sumLength,[1,numRed]);


%% Test
memOther = zeros(numRep,maxSyn+1);
memRed = zeros(numRep,maxSyn+1);
memRat = zeros(numRep,11);
testStat = zeros(numRep,1); %test Statistic 1

for r = 1:numRep
    
    randOther = rand(1,numOther);
    randRed = rand(1,numRed);
    
    matOther = repmat(randOther,[axNum 1]);
    matRed = repmat(randRed, [axNum 1]);
    
    testOther = otherSumLength<matOther;
    testRed = redSumLength<matRed;
    
    resultOther  = sum(testOther,1)+1;
    resultRed = sum(testRed,1)+1;
    
    histOtherAx = hist(resultOther,1:1:axNum);
    histRedAx = hist(resultRed,1:1:axNum);
    
    scatter(axLength,histRedAx,'.'),pause(.1)
    
    
    
    ratRedOther = histRedAx./(histOtherAx+histRedAx);
    ratRedOther(isnan(ratRedOther)) = 0;
    
    histOther = hist(histOtherAx,0:1:maxSyn);
    histRed = hist(histRedAx,0:1:maxSyn);
    histRat = hist(ratRedOther,0:.1:1);
    
    
    memOther(r,:) = histOther;
    memRed(r,:) = histRed;
    memRat(r,:) = histRat;
    
    
    %%Get test statistic
    distFromAv = ratRedOther - averageRat;
    testStat = mean(distFromAv(distFromAv>0));
   
    
end
    
meanOther = mean(memOther,1);
meanRed = mean(memRed,1);
meanRat = mean(memRat,1);

%% Display

subplot(4,2,2)
bar([0:1:maxSyn],meanOther)
title('Test Other')
subplot(4,2,4)
bar([0:1:maxSyn],meanRed)
title('Test Red')
subplot(4,2,6)
bar([0:.1:1],meanRat)
title('Test Ratio')


%% Analyze Results
P = sum(sumDistFromAv>=realSumDistFromAv)/numRep;
disp(sprintf('P = %f',P))

subplot(4,2,[7 8])
hist(sumDistFromAv);
hold on
scatter(realSumDistFromAv,1,'r')
xlim([0 1])
hold off
    
    
    
    

