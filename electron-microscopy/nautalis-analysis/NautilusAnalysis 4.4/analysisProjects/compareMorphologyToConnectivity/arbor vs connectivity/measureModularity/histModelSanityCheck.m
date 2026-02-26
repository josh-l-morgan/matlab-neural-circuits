
clear all
MPN = GetMyDir
fileName = sprintf('%sskelOverlapPred.mat',MPN);
load(fileName,'skelOverlapPred')
seedList = [108 201];

axCellHists = skelOverlapPred.axCellHists;
histBin = skelOverlapPred.histBin;
conDists = skelOverlapPred.conDists;
sumHists = axCellHists{1,1} * 0;
allCons = [];
for y = 1:size(axCellHists,1)
    for x = 1:size(axCellHists,2)
        
        sumHists = sumHists + axCellHists{y,x};
        allCons = cat(1,allCons,conDists{y,x});
    end
end
bar(histBin,sumHists)
conHist = hist(allCons,histBin);
maxCon = max(find(conHist>0));
maxConDist = histBin(maxCon)

hold on
bar(histBin,conHist,'r')
hold off
