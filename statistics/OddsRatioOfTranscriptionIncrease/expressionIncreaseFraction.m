

%% Set variables
ng1 = 20000; % number of genes measured in dataset 1
ng2 = 20000; % number of genes measured in dataset 2
pickTop = 100; % number of top genes selected from dataset 1
reps = 100000; % number of times to run simulation
CI = 0.95;

%% Make test data
realDat1 = randn(ng1,1); % Make dataset 1;
[sort1 idx] = sort(realDat1,'descend'); %Sort results
realPick = idx(1:pickTop); %Find top pickTop results

%% Make data set 2
realDat2 = randn(ng2,1); % random values;
realDat2 = mean([realDat1 realDat2],2); %mix with dataset 1

%% Measure
realVals1 = realDat1(realPick); %Get top values from group 1
realVals2 = realDat2(realPick); %Get top values from group 2
realPosFrac1 = mean(realVals1>0); %Measure number of top values above 0 
realPosFrac2 = mean(realVals2>0); %Measure number of top values above 0 


randPosFrac = zeros(reps,1);
for r = 1:reps % repeat simulation

    randPick = randsample(ng2,pickTop);
    randVals2 = realDat2(randPick);
    randPosFrac(r) = mean(randVals2>0);

end


%% Calculate
sortRand = sort(randPosFrac,'ascend');
lowBound = sortRand(round(r*(1-CI)/2));
highBound = sortRand(round(r*(1-(1-CI)/2)));

resStr = sprintf('%0.1f%% increased (%0.1f%% to %.01f%% in random selections)',...
    realPosFrac2*100,lowBound*100, highBound*100)

%% Display
clf
hRange = [0:.01:1];
histRandPosFrac = hist(randPosFrac,hRange)/r;
maxY = max(histRandPosFrac(:));
plot(hRange,histRandPosFrac,'k')
hold on
plot([lowBound lowBound],[0 maxY],'r')
plot([highBound highBound],[0 maxY],'r')
scatter(realPosFrac2,0,'r','filled')
title(resStr)
xlabel('fraction that increased')













