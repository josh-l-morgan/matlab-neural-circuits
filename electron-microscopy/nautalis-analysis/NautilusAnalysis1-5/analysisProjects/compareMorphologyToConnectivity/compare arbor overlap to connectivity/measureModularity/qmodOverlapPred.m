
clear all
load('MPN.mat')
fileName = sprintf('%sskelOverlapPred.mat',MPN);
load(fileName,'skelOverlapPred')
seedList = [108 201];

useList = skelOverlapPred.useList;
useList.seedList = seedList;
Qreal = use2Modularity(useList);
[shuffleGraph,newOrder,Q]  = SortPlusUse2Modularity(useList);


skelOverlapPred.axCellHists

conRaw = useList.con;
mask = conRaw*0+1;
if 1 %zero seeds
    for i = 1:length(seedList)
        mask(:,find(useList.postList == seedList(i))) = 0;
    end
end

binRes = .1;
[predSyn nodeFrac] = overlapPredictSynOfSubset(skelOverlapPred,binRes,mask);
%con = round(predSyn);
%con = conRaw;
[predCon] = predNewCon(conRaw,predSyn,mask);
predCon(~mask) = conRaw(~mask);

%%
useRange = useList;

resRange = [0.1: 1:100];
QrangeSort = zeros(0,length(resRange));
QrangeFrac = QrangeSort;
QrangeRound = QrangeSort;
for i = 1:length(resRange)
    [predSyn nodeFrac] = overlapPredictSynOfSubset(skelOverlapPred,resRange(i),mask);
    %con = round(predSyn);
    %con = conRaw;
    [predCon] = predNewCon(conRaw,predSyn,mask);
    predCon(~mask) = conRaw(~mask);
    useRange.con = predCon;
    QrangeSort(i) = use2Modularity(useRange);
    
    
    [predSyn synNum] = overlapFracPredictSynOfSubset(skelOverlapPred,resRange(i),mask);
    
    
    [predCon] = predNewCon(conRaw,predSyn,mask);
    predCon(~mask) = conRaw(~mask);
    useRange.con = predCon;
    QrangeFrac(i) = use2Modularity(useRange);
    
    predCon = round(predSyn);
    predCon(~mask) = conRaw(~mask);
    useRange.con = predCon;
    QrangeRound(i) = use2Modularity(useRange);
    
end


[Q, shuffleGraph,newOrder]  = SortPlusUse2Modularity(useRange);

plot(resRange,QrangeSort,'b')
hold on
plot(resRange,QrangeFrac,'r')
plot(resRange,QrangeRound,'g')

hold off


useRes = (resRange > 8) & (resRange <= 100);
predVals = QrangeFrac(useRes);
mean(predVals)
min(predVals)
max(predVals)


%%
binRes = .1;
[predSyn nodeFrac] = overlapPredictSynOfSubset(skelOverlapPred,binRes,mask);
%con = round(predSyn);
%con = conRaw;
[randCon] = randNewCon(conRaw,rand(size(conRaw)),mask)
randCon(~mask) = conRaw(~mask);




usePred = useList;
usePred.con = predCon;
useRand = useList;
useRand.con = randCon;


Qreal = use2Modularity(useList)
Qpred = use2Modularity(usePred)
Qrand = use2Modularity(useRand)




