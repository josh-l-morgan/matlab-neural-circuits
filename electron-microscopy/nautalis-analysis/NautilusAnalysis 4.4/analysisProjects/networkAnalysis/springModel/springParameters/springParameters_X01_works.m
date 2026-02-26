function[springDat] = springParameters(springIn)


%% Create springDat (need nodeCol, nodeIDs, 

nodeIDs = springIn.nodeIDs;
allEdges = springIn.allEdges;
nodeCol = springIn.nodeCol;
nodeType = springIn.nodeType;
seedList = springIn.seedList;
allWeights = springIn.allWeights;

nodes.type = nodeType;

seed.list = seedList;
% seed.seedPref = synPref;
% seed.isPref = isPref;


%% Get new edges
[edges] = uniqueEdges(allEdges,nodeIDs); 

edgeCol = nodeCol(edges.e1,:);
edgeCol(edges.ew==0,:) = edgeCol(edges.ew==0,:)*.2;

edgeCol(edgeCol<0) = 0;
edgeCol(edgeCol>1) = 1;

edges.col = edgeCol;
edges.symmetry = edges.ew*0;
edges.width = sqrt(edges.ew);
edges.text = 0;
%edges.width = edges.ew*1;


nodes.labelIDs = nodeIDs;
nodes.color = nodeCol;

%% set param
param.k = .1;
param.disperse = 1;
param.centerSpring = 1;


param.reps = 4000;
param.zeroSeeds = 1;
param.noise = (param.reps:-1:1)/param.reps*1;
param.damp = .5;
param.fsize = 200;
param.speedMax = 1;
param.imageFreq = 100;

param.startRepulse = param.reps -500;%
param.repulse = 10; %1
param.minDist = 5;

%%Make groups
%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');

groupNum = max(nodeType);

group(1).ind = find(nodeType==1);
group(2).ind = find(nodeType==2);
group(3).ind = find(nodeType==3);
group(4).ind = find(nodeType==4);

group(1).marker = '^';
group(2).marker = 'o';
group(3).marker = '^';
group(4).marker = 'o';

group(1).size = 20;
group(2).size = 20;
group(3).size = 150;
group(4).size = 150;

group(1).lineWidth = .5;
group(2).lineWidth = .5;
group(3).lineWidth = .5;
group(4).lineWidth = .5;

group(1).text = 0;
group(2).text = 0;
group(3).text = 0;
group(4).text = 0;


springDat.group = group;
springDat.groupNum = groupNum;
springDat.param = param;
springDat.nodes = nodes;
springDat.edges = edges;
springDat.seed = seed;

