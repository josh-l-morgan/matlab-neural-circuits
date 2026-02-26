function[springDat] = skelSpringParameters_01(nep)



%% translate for spring dat
edges = uniqueEdgesSimple(nep.edges,nep.nodes); 
%edgeCol(edges.ew==0,:) = edgeCol(edges.ew==0,:)*.2;
edges.col = nep.edgeCol;
edges.symmetry = edges.ew*0;
edges.width = nep.edgeWidth; %[.2*(edges.ew).^(1/2)];
edges.text = 0;

nodes.labelIDs = nep.nodes;
nodes.color = nep.nodeCol;
nodes.mass = nep.nodeMass;
nodes.pushes = nep.nodePushes;
nodes.isPushed = nep.nodeIsPushed;
nodes.type = nep.nodeType;



%% set param
param.k = 1;
param.disperse = 1;
param.centerSpring = 1;


param.reps = 10;
param.zeroSeeds = 1;
param.noise = (param.reps:-1:1)/param.reps*1;
param.damp = .5;
param.fsize = 200;
param.speedMax = 1;
param.imageFreq =100;

param.startRepulse = param.reps -500;%
param.repulse = 10; %1
param.minDist = 5;

%% Make groups
%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');

groupNum = max(nep.nodeType);

group(1).ind = find(nep.nodeType==1);
group(2).ind = find(nep.nodeType==2);
group(3).ind = find(nep.nodeType==3);
group(4).ind = find(nep.nodeType==4);

group(1).marker = 'o';
group(2).marker = 'o';
group(3).marker = 'o';
group(4).marker = 'o';

group(1).size = 30;
group(2).size = 2;
group(3).size = 2;
group(4).size = 300;

group(1).lineWidth = .1;
group(2).lineWidth = .1;
group(3).lineWidth = .1;
group(4).lineWidth = .1;

group(1).text = 0;
group(2).text = 0;
group(3).text = 0;
group(4).text = 0;

%% make springDat
springDat.group = group;
springDat.groupNum = groupNum;
springDat.param = param;
springDat.nodes = nodes;
springDat.edges = edges;
springDat.seed.list = nep.seedList;
springDat.nep = nep;

