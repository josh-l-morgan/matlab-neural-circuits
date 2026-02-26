
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])
allEdges = obI.nameProps.edges(:,[2 1]);

useList = obI2cellList_seedInput_RGC_TCR(obI,[108 201 903 907]);
useList = obI2cellList_seedInput_RGC_TCR(obI,[108 201]);

targ = find(useList.preList == 1121);
notTarg = setdiff([1:length(useList.preList)],targ);
useList.preList = useList.preList(notTarg);
useList.con = useList.con(notTarg,:);

conPref = seedPreferences([108 201],useList);



crossoverAxons = [2032	2033	2034	2035]

%%
clf

cellNum = length(conPref.cellList);

con = useList.con
targ = find(useList.postList == 108);
isPre = con(:,targ)>0;
linkedSyn = con(isPre,:)
sum(linkedSyn(:))-sum(con(:,targ))
sum(sum(linkedSyn,1)>0)

%%
clf

minAx = min(conPref.sharedAx,[],1);
sumAx = sum(conPref.sharedAx,1);
linkedSum = sum(conPref.sharedAx>0,2)-1;
multSums = sumAx(minAx>0);
countCrossed = sum(minAx>0);
histBin = [1:max(sumAx)];
histMin = hist(multSums,histBin);
histSum = hist(sumAx,histBin);
bar(histBin,histSum,'b')
hold on
bar(histBin,histMin,'r')
hold off


multNonSeed = sum(sumAx>1)-2;
crossoverNum = sum(minAx>0);
PercentCrossoverOfMults = crossoverNum/multNonSeed * 100

PercentCrossover = sum(minAx>0)/length(minAx) * 100;
disp(sprintf('Percent Crossover %.1f%%',PercentCrossover))


scatter(conPref.sharedSyn(1,:),conPref.sharedSyn(2,:))


%%

is108 = conPref.sharedAx(1,:)>0

is108 & (sumAx>1)


%%  Compare "crossed to uncrossed"

isMixed = minAx>0;
pre201 = preTo(allEdges,201);
unCrossed = intersect(pre201(:,1),useList.preList);
unCrossed = setdiff(unCrossed,crossoverAxons);

unXmix = zeros(length(unCrossed),2);
for i = 1:length(unCrossed)
    postI = postTo(allEdges,unCrossed(i));
    mixCount = zeros(size(postI,1),1);
    for p = 1:size(postI,1)
        targ = find(useList.postList == postI(p));
        if ~isempty(targ)
            unXmix(i,1) = unXmix(i,1) + isMixed(targ);
            unXmix(i,2) = unXmix(i,2) + 1;
        end
    end
end

Xmix = zeros(length(crossoverAxons),2);
for i = 1:length(crossoverAxons)
    postI = postTo(allEdges,crossoverAxons(i));
    mixCount = zeros(size(postI,1),1);
    for p = 1:size(postI,1)
        targ = find(useList.postList == postI(p));
        if ~isempty(targ)
            Xmix(i,1) = Xmix(i,1) + isMixed(targ);
            Xmix(i,2) = Xmix(i,2) + 1;
        end
    end
end

unXrat = unXmix(:,1)./unXmix(:,2);
Xrat = Xmix(:,1)./Xmix(:,2)
ranksum(unXmix(:,1)./unXmix(:,2) , Xmix(:,1)./Xmix(:,2))

sprintf('%.3f',[unCrossed' unXmix unXrat])
[unCrossed' unXmix unXrat*100]
[crossoverAxons' Xmix Xrat*100]


median(unXrat)
median(Xrat)
rangeX(unXrat)
rangeX(Xrat)
ranksum(Xrat,unXrat)

%% cellA

isMixed = minAx>0;
pre108 = preTo(allEdges,108);
cont108 = intersect(pre108(:,1),useList.preList);

contMix = zeros(length(cont108),2);
for i = 1:length(cont108)
    postI = postTo(allEdges,cont108(i));
    mixCount = zeros(size(postI,1),1);
    for p = 1:size(postI,1)
        targ = find(useList.postList == postI(p));
        if ~isempty(targ)
            contMix(i,1) = contMix(i,1) + isMixed(targ);
            contMix(i,2) = contMix(i,2) + 1;
        end
    end
end

contRat = contMix(:,1)./contMix(:,2);
[cont108' contMix contRat*100]





%%
postUnMix = [];
for i = 1:length(unCrossed)
    postI = postTo(allEdges,unCrossed(i));
    postUnMix = [postUnMix postI(:,1)'];
end
postUnMix = unique(postUnMix);
postUnMix = intersect(useList.postList,postUnMix);

allUnXmix = [];
for p = 1:length(postUnMix)
        targ = find(useList.postList == postUnMix(p));
        if ~isempty(targ)
            allUnXmix = [allUnXmix  isMixed(targ)];
        end
end

postMix = [];
for i = 1:length(crossoverAxons)
    postI = postTo(allEdges,crossoverAxons(i));
    postMix = [postMix postI(:,1)'];
end
postMix = unique(postMix);
postMix = intersect(useList.postList,postMix);

allXmix = [];
for p = 1:length(postMix)
        targ = find(useList.postList == postMix(p));
        if ~isempty(targ)
            allXmix = [allXmix  isMixed(targ)];
        end
end

mean(allXmix)
mean(allUnXmix)
ranksum(allUnXmix,allXmix)
    
    
    
    
    
    
    
    
    
    
    
    

