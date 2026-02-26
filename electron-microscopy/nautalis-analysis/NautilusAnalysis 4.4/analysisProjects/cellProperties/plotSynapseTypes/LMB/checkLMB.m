%% Load data

clear all
load('MPN.mat')
%MPN = GetMyDir;
load([MPN 'obI.mat']);

seedList = [ 108  201 907 903];

useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
seedPref = seedPreferences(seedList,useList);
allEdges = obI.nameProps.edges(:,[2 1]);


%% Find ldm
sum(obI.nameProps.ldm)
sum(obI.nameProps.LMB)

obI.nameProps.names{obI.nameProps.ldm}
[obI.nameProps.allIDs{obI.nameProps.ldm}]






































