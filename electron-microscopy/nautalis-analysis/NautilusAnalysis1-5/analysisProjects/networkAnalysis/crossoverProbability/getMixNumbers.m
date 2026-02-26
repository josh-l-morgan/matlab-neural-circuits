
%%Test the probability that crossover synapses would be concentrated on few
%%axons
clear all
load('MPN.mat')
load([MPN 'obI.mat'])
allEdges = obI.nameProps.edges(:,[2 1]);
seedList = [108 201 903 907];
useCells = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);

edges = unique(allEdges,'rows');


useTCR = useCells.postList;


isA = (intersect(useTCR,conTo(1).tcrList))
isB = (intersect(useTCR,conTo(2).tcrList))
isD = (intersect(useTCR,conTo(4).tcrList))

isAD = intersect(isA,isD)
isBD = intersect(isB,isD)



