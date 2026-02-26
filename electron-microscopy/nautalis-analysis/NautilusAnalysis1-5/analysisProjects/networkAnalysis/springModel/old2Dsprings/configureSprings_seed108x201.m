
% %%
% F = -kX
% k = stiffness
% X = proportional to distance;

%% Get data
loadData = 1;
if loadData
    clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    seedList = [108 201];
    useList = obI2cellList_seedInput(obI,seedList);
    %useList = obI2cellList_all(obI);
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
end

springDir = 'D:\LGNs1\Analysis\springDat3\'
if ~exist(springDir,'dir'), mkdir(springDir), end


[someEdges ew] = uniqueEdges(allEdges,nodeIDs);


%% shape attributes


postSynPref = seedPref.sharedSynNorm(1,:)./sum(seedPref.sharedSynNorm,1);
preSynPref = seedPref.ax2seed(1,:)./sum(seedPref.ax2seed,1);
synPref = [preSynPref postSynPref];
%synPref = preSynPref;
grouping = [preSynPref * 0 + 1   postSynPref * 0 + 2];
isPref = ~isnan(synPref);

nodeNum = length(isPref);
nodeCol = zeros(nodeNum,3);
nodeCol(isPref,1) = synPref(isPref);
nodeCol(isPref,3) = 1-synPref(isPref);
nodeCol(~isPref,2) = 1;
%set(0,'DefaultAxesColorOrder',nodeCol)

seed.list = seedList;
seed.seedPref = synPref;
seed.isPref = isPref;


%% edge data
edges.all = allEdges;

rawCon = zeros(nodeNum,nodeNum);
useCon = rawCon * 0;
for i = 1:length(nodeIDs)
    for p = 1:length(nodeIDs)
        rawCon(i,p) = rawCon(i,p) + sum((allEdges(:,1) == nodeIDs(i)) & (allEdges(:,2) == nodeIDs(p)));
        rawCon(i,p) = rawCon(i,p) + sum((allEdges(:,1) == nodeIDs(p)) & (allEdges(:,2) == nodeIDs(i)));
        useCon(i,p) = (p>i);
    end
end
con  = rawCon;

%%Get new edges
[e1 e2] = find((con.*useCon)>0);
ei = sub2ind(size(con),e1,e2);
ew = con(ei);
edgeCol = nodeCol(e1,:) + nodeCol(e2,:);
edgeCol(ew==0,:) = edgeCol(ew==0,:)*.2;
edgeCol(edgeCol<0) = 0;
edgeCol(edgeCol>1) = 1;

edges.all = allEdges;
edges.e1 = e1;
edges.e2 = e2;
edges.ew = ew;
edges.col = edgeCol;
edges.symmetry = ew*0+1;
%edges.width = edgeWidth;




%% Node data

nodes.labelIDs = [useList.preList useList.postList];
nodes.type = zeros(nodeNum,1);
%     nodes.grouping1(nodeIDs<1000) = 1;
%     nodes.grouping1((nodeIDs>=1000) & (nodeIDs < 10000)) = 2;
%     nodes.grouping1(nodeIDs>=10000) = 3;
%
nodes.grouping1 = grouping;
nodes.color = nodeCol;


%%set param
param.zeroSeeds = 1;
param.k = 1;
param.noise = [10 1 .1 0];
param.damp = .5;
param.disperse = 5;
param.reps = 2000;
param.fsize = 200;
param.speedMax = 1;
param.centerSpring = 2;
param.imageFreq = 1;

param.startRepulse = 500;%
param.repulse = 5; %1
param.minDist = 5;

%%Make groups
%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');

group(1).ind = nodes.grouping1 == 1;
group(2).ind = nodes.grouping1 == 2;
group(3).ind = nodes.grouping1 == 3;
group(1).marker = '^';
group(2).marker = 'o';
group(3).marker = 's';
group(1).size = 100;
grooup(2).size = 5;
group(3).size = 20;
groupNum = length(group);


springDat.group = group;
springDat.param = param;
springDat.nodes = nodes;
springDat.edges = edges;
springDat.seed = seed;



results = runSprings(springDat);

%results = runSpringMass(springDat);







