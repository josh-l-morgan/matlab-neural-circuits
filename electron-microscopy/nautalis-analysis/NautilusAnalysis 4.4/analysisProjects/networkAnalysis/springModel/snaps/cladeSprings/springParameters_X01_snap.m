function[springDat] = springParameters(springIn)


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
edges.width = [sqrt(syns.ew); juncs.ew * 3];

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
param.disperse = .01; %default = 1
param.centerSpring = 1;


param.reps = 60000;
repCount = 1:param.reps;
param.zeroSeeds = 0;
param.noise = (repCount<1000)*1 + (repCount>1000) * .01;%(param.reps:-1:1)/param.reps*1 + .01;
param.damp = .5;
param.fsize = 200;
param.speedMax = 1;
param.imageFreq = 100;

param.startRepulse = param.reps ;%
param.repulse = 0; %1
param.minDist = 0;

param.doSnap = 1;
param.snap.start = 3000;
param.snap.wait = 100;
param.snap.longest = 1;
param.snap.synAll = 0;

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

