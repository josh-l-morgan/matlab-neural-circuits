


%% Plot axon to tcr
clear all
MPN = GetMyDir;
load([MPN 'obI.mat'])
shouldPause = 10;
anchorScale = [-.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8 * 4 anchorScale(2) * 8 * 4 anchorScale(3)* 4 * 4];


%% variables

shiftList = [-20 -10 -5 -4 -3 -2 -1.5 -1 -.5 -.25 0 .25 .5 1 1.5 2 3 4 5 10 20]

skelOverlapPred.predictionType = 'skelOverlap';
skelOverlapPred.shiftList = shiftList;

firstPredThresh = 0;
binWidth = .1;
skelOverlapPred.binWidth = binWidth;
seedList = [108 201]
crossoverAxons = [2032	2033	2034	2035]
noSkel = [2014 1026]

% muchTraced = [106 107 108 109 111 112 117 120 123 129 133 134 148 156 159 162 163 169 ...
%     170 201 203 205 206 207 210 212 213 215 216 218];
% skelOverlapPred.muchTraced = muchTraced;



%% Get cells

useList = obI2cellList_seedInput(obI,seedList);
axList = useList.preList;
cellList = useList.postList;
synMat = useList.con;
axList = setdiff(axList,noSkel);


useList.preList = axList;
skelOverlapPred.useList = useList;

synapses = obI.nameProps.edges;
edges = synapses(:,1:2);


disp(sprintf('Results calculated without %d',noSkel));

%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end




%%  Show 108 and non-seed 108 synapse histograms
targ = find(cellList == 108)
seedSyn = con(:,targ);
pre108 = find(seedSyn>0);

for i = 1:length(pre108)
   linkedToPre = setdiff(find(con(pre108(i),:)>0),targ);
   synNums{i}  = con(i,linkedToPre);
end

allSynNums = [synNums{:}];

synBin = 1:15;
histSeed = hist(con(pre108,targ),synBin);
histNonSeed = hist(allSynNums,synBin)

bar([histSeed; histNonSeed]')

