 

clf
maxDist = .025;
medPred = zeros(length(POI.testDepths),length(testLengths));

d  = 16

bestPred = usePredN(:,d,round(bs1m),round(bs2m),round(bnm));

polRange = [-1:.05:1];
for p = 1:length(POI.testDepths)

    subplot(length(POI.testDepths),1,p)
    isPlane = (POI.plane(useRoi) == p) & (POI.planeDist(useRoi)<=maxDist);

    pols = bestPred(isPlane);
    medPred(p,d) = median(pols);
    histPols = hist(pols,polRange);
    bar(polRange,histPols,'facecolor','k','edgecolor','none')
%     title(sprintf('length constant %0.2f, testDepth = %0.2f',...
%         testLengths(d),POI.testDepths(p)))
    title(sprintf('testDepth = %0.2f', POI.testDepths(p)))
    ylim([0 20])
end

