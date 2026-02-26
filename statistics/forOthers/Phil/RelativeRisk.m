

% dat1 = []
% dat2 = []

rats = cat(1,dat1,dat2);
dies = cat(1,dat1*0,dat2*0+1);

%% reality check
testThresh = 1.3;
isActive = rats>=testThresh;
activeLives = sum(dies(isActive)==0);
activeDies = sum(dies(isActive)==1);
activeSurvivalRate = activeLives/(activeLives+activeDies)

isNotActive = rats<testThresh;
notActiveLives = sum(dies(isNotActive)==0);
notActiveDies = sum(dies(isNotActive)==1);
notActiveSurvivalRate = notActiveLives/(notActiveLives+notActiveDies)


%%

threshes = [min(rats):.01:max(rats)];
groupSizeThresh = 0.05; %maller group must be X fraction of total

totNum = length(rats);
minSize = ceil(totNum * groupSizeThresh);

clear numIn numOut numInDies numOutDies 
for t = 1:length(threshes)
   
    isIn = find(rats>=threshes(t));
    isOut = find(rats<threshes(t));
    
    numIn(t) = length(isIn);
    numOut(t) = length(isOut);
    numInDies(t) = sum(dies(isIn));
    numOutDies(t) = sum(dies(isOut));
    
    
end

inRisk = numInDies./numIn;
inOdds = numInDies./(numIn-numInDies);
outRisk = numOutDies./numOut;
outOdds = numOutDies./(numOut-numOutDies);

smaller = min(numIn, numOut);
isValid = find(smaller>0);
isValid2 = find(smaller>=minSize);
validBound = [threshes(isValid2(1)) threshes(isValid2(end))];

relativeRisk = inRisk./outRisk;
oddsRatio = inOdds./outOdds;

subplot(4,1,1)

plot(threshes(isValid),numIn(isValid),'r')
hold on
title('Number in high calcium (red) and low calcium(green) at each threshold')
plot(threshes(isValid),numOut(isValid),'g')
hold off

subplot(4,1,2)
plot(threshes(isValid),inRisk(isValid),'r')
hold on
title('risk (number that die / total number in group')
plot(threshes(isValid),outRisk(isValid),'g')
hold off

subplot(4,1,3)
plot(threshes(isValid),relativeRisk(isValid),'k')
hold on
title('relative risk.  Risk of above calcium threshold group/ risk of others')
plot(threshes(isValid),oddsRatio(isValid),'m')
plot([validBound(1) validBound(1)],[0 max(relativeRisk(isValid))],'r')
plot([validBound(2) validBound(2)],[0 max(relativeRisk(isValid))],'r')
hold off

subplot(4,1,4)
plot(threshes(isValid2),relativeRisk(isValid2),'k')
hold on
title('relative risk for ranges in which smalles group size is big enough')
plot(threshes(isValid2),oddsRatio(isValid2),'m')
hold off


