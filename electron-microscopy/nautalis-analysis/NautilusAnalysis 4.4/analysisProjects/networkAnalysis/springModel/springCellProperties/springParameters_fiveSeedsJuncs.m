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


%%Get new edges
[syns] = uniqueEdges(allEdges,nodeIDs); 
juncs = uniqueEdges(sort(allJuncs','descend')',nodeIDs);

edgeCol = nodeCol(syns.e1,:);
edgeCol(syns.ew==0,:) = edgeCol(syns.ew==0,:)*.2;

edgeCol(edgeCol<0) = 0;
edgeCol(edgeCol>1) = 1;

edges.col = cat(1,edgeCol,ones(length(juncs.ew),3));
edges.symmetry = [syns.ew; juncs.ew]*0;
edges.width = [sqrt(syns.ew); juncs.ew * 1];

edges.e1 = [syns.e1;juncs.e1];
edges.e2 = [syns.e2;juncs.e2];
edges.ew = [syns.ew;juncs.ew * 0];
edges.use = [syns.use;juncs.use];
edges.all = [syns.all;juncs.all];





nodes.labelIDs = nodeIDs;
nodes.color = nodeCol;

%%set param
param.k = .1;
param.disperse = 1.5;
param.centerSpring = 1;


param.reps = 4000;
param.zeroSeeds = 0;
param.noise = (param.reps:-1:1)/param.reps*1;
param.damp = .5;
param.fsize = 200;
param.speedMax = 1;
param.imageFreq = 100;

param.startRepulse = param.reps -500;%
param.repulse = 10; %1
param.minDist = 6;

%%Make groups
%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');

groupNum = 3;

group(1).ind = find(nodeType==1);
group(2).ind = find(nodeType==2);
group(3).ind = find(nodeType==3);
group(1).marker = '^';
group(2).marker = 'o';
group(3).marker = 'o';
group(1).size = 100;
group(2).size = 50;
group(3).size = 20;
group(1).lineWidth = .5;
group(2).lineWidth = .7;
group(3).lineWidth = .7;

springDat.group = group;
springDat.groupNum = groupNum;
springDat.param = param;
springDat.nodes = nodes;
springDat.edges = edges;
springDat.seed = seed;

