function[springDat] = springParameters_X01_figDraftNoProp(springIn)


%% Create springDat (need nodeCol, nodeIDs, 

nodeIDs = springIn.nodeIDs;
allEdges = springIn.allEdges;
nodeCol = springIn.nodeCol;
nodeType = springIn.nodeType;
seedList = springIn.seedList;
allWeights = springIn.allWeights;

try allJuncs = springIn.allJuncs;
catch err, allJuncs = []; end

nodes.type = nodeType;

seed.list = seedList;
% seed.seedPref = synPref;
% seed.isPref = isPref;


%% Get new edges
syns = uniqueEdges(allEdges,nodeIDs); 
juncs = uniqueEdges(sort(allJuncs','descend')',nodeIDs);

edgeCol = nodeCol(syns.e1,:);
edgeCol(syns.ew==0,:) = edgeCol(syns.ew==0,:)*.2;


edges.col = cat(1,edgeCol,ones(length(juncs.ew),3)*1000);
edges.symmetry = [syns.ew; juncs.ew]*0;
edges.width = [2*(syns.ew).^(1/2); juncs.ew * 3];

edges.e1 = [syns.e1;juncs.e1];
edges.e2 = [syns.e2;juncs.e2];
edges.ew = [syns.ew;juncs.ew * 0];
edges.use = [syns.use;juncs.use];
edges.all = [syns.all;juncs.all];

edges.text = 0;
%edges.width = edges.ew*1;


nodes.labelIDs = nodeIDs;
nodes.color = nodeCol;

%% set param
param.k = .1;
param.disperse = 10;
param.centerSpring = 1;


param.reps = 4000;
param.zeroSeeds = 1;
param.noise = (param.reps:-1:1)/param.reps*1;
param.damp = .5;
param.fsize = 200;
param.speedMax = 1;
param.imageFreq = 10;

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

group(1).size = 80;
group(2).size = 80;
group(3).size = 600;
group(4).size = 600;

group(1).lineWidth = 1;
group(2).lineWidth = 1;
group(3).lineWidth = 1;
group(4).lineWidth = 1;

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

