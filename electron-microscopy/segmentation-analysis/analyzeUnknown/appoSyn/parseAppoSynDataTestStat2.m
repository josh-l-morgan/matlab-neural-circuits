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

probIdx = sub2ind(size(dat),realCountAp,realCountAx+1);
realCumProbs = cumProb(probIdx);
realMeanProb = mean(realCumProbs);

realMultSyn = sum(realCountAx-1);

realHistAx = hist(realCountAx,histBin)


%%  Model
reps = 10000;
modSyn = realSyn*0;
apNum = length(ids);
synNum = sum(realSyn);

histAx = zeros(reps,6);
for r = 1:reps
    modSyn = modSyn*0;
    
    genRand = rand(1,apNum);
    [sortRand idx] = sort(genRand);
    pickRand = idx(1:synNum);
    modSyn(pickRand) = 1;
    trackSyn(r) = sum(modSyn);
    
    %modSyn = realSyn;
    
    withSyn = ids(modSyn>0);
    uAx = unique(withSyn);
    clear countAx countAp
    for i = 1:length(uAx)
        countAx(i) = sum(withSyn==uAx(i));
            countAp(i) = id2ap(uAx(i));
    end
    
    
probIdx = sub2ind(size(cumProb),countAp,countAx+1);
cumProbs = cumProb(probIdx);
meanProb(r) = mean(cumProbs);
    
    multSyn(r) = sum(countAx-1);
    trackAx(r) = sum(countAx);
    
    histAx(r,:) = hist(countAx,histBin);
    
    if ~mod(r,1000)
        disp(sprintf('running %d of %d',r,reps))
    end
    
end





%%
subplot(1,1,1)
histHistBin = [1.1:80.1];
clear lowC midC highC Xval realVal
for i = 1:3
    h = i+1;
    
    
    hDat = histAx(:,h);
    sortDat = sort(hDat,'ascend');
    lowC(i) = sortDat(round(.025*reps));
    midC(i) = sortDat(round(.5*reps));
    highC(i) = sortDat(round(.975*reps));
    meanC(i) = mean(hDat);
    
    Xval(i) = i;
    realVal(i) = realHistAx(h);
    
    
    
    hist3 = hist(histAx(:,h),histHistBin);
    subplot(3,2,(i-1)*2+1)
    bar(histHistBin,hist3/reps,'r')
    hold on
    scatter(realHistAx(h),0,40,'b','filled')
    title(sprintf('Axons with %d synapses',i))
    hold off
    
    
    difC(:,i) = histAx(:,h)-midC(i);
    realDifC(:,i) = realHistAx(:,h)-midC(i);
    
    
end

subplot(3,2,[2 4 6])

bar(Xval,midC,'FaceColor','none','EdgeColor','r','BarWidth',.9)
hold on

bar(Xval,realVal,'b','BarWidth',.6)

errorbar(Xval,midC,midC-lowC,midC-highC,'Color','k','lineWidth',2,...
    'LineStyle','none')
scatter(Xval,midC,40,'r','filled')
hold off


%% calc P for all

sumDifs = sum(abs(difC),2);
realSumDif = sum(abs(realDifC),2);

P = sum(sumDifs>=realSumDif)/reps;

sortDifs = sort(sumDifs,'ascend');
dif95 = sortDifs(round(.95*reps));
difMid = median(sumDifs);

text(2,60,sprintf('normal dif = %.1f', difMid))
text(2,50,sprintf('95 dif = %.1f',dif95))
text(2,40,sprintf('real dif = %.1f',realSumDif))
text(2,30,sprintf('P = %.7f',P))
title('Model (red) vs Observed (blue)')



%% calc P for 1
% 
% sumDifs = abs(difC(:,2));
% realSumDif = abs(realDifC(2));

