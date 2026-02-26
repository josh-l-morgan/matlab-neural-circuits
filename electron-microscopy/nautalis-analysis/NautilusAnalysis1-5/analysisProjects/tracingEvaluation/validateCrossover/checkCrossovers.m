

%% Get crossover TCR / RGC
targCells = [108 201];
conTo = makeConTo(obI,targCells);
allEdges = obI.nameProps.edges(:,[2 1]);



crossTCR = intersect(conTo(1).tcrList, conTo(2).tcrList)

for i = 1:length(crossTCR)
    pre2 = preTo(allEdges,crossTCR(i));
    preAx = pre2(:,1)';
    crossAx{i,1} = intersect(conTo(1).rgcList,preAx);
    crossAx{i,2} = intersect(conTo(2).rgcList,preAx);
end






























