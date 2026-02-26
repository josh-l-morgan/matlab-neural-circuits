function[springDat] = skelSpringParameters_01(nep)



%% translate for spring dat

conNodes = unique(nep.edges(:));
conNodes = intersect(nep.nodes,conNodes);

useNodes = zeros(length(nep.nodes),1)>0;
useNodes(conNodes) = 1;

nodes.labelIDs = nep.nodes(useNodes);
nodes.nodeName = nep.nodeName(useNodes);
nodes.color = nep.nodeCol(useNodes,:);
nodes.mass = nep.nodeMass(useNodes);
nodes.pushes = nep.nodePushes(useNodes);
nodes.isPushed = nep.nodeIsPushed(useNodes);
nodes.type = nep.nodeType(useNodes);
nodes.usePos = nep.usePos(useNodes);
nodes.nodePos = nep.nodePos(useNodes,:);
nodes.nodeMove = nep.nodeMove(useNodes);

nodes.group = nodes.type;

%edges = uniqueEdgesSimple(nep.edges,conNodes); 
[useEdge newEdges] = usedEdges(nep.edges,conNodes); 

%edgeCol(edges.ew==0,:) = edgeCol(edges.ew==0,:)*.2;
edges.use = newEdges;
edges.all = nep.edges(useEdge,:);
edges.e1 = newEdges(:,1);
edges.e2 = newEdges(:,2);
edges.ew = nep.edgeWeights(useEdge);
edges.col = nep.edgeCol(useEdge,:);
edges.symmetry = edges.ew*0;
edges.width = nep.edgeWidth(useEdge); %[.2*(edges.ew).^(1/2)];
edges.text = 0;





%% set param
param.k = 1;
param.disperse = 1;
param.centerSpring = 0.01;


param.reps = 10000;
param.usePos = 1;
param.zeroSeeds = 1;
param.noise = (param.reps:-1:1)/param.reps*1;
param.damp = .5;
param.fsize = 200;
param.speedMax = 1;
param.imageFreq =100;


param.startRepulse = param.reps -500;%
param.repulse = 10; %1
param.minDist = 5;

param.edgeColor = [.5 .5 .5]    ;
%% Make groups
%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');

groupNum = max(nep.nodeType);

group(1).ind = find(nodes.type==1);
group(2).ind = find(nodes.type==2);
group(3).ind = find(nodes.type==3);
group(4).ind = find(nodes.type==4);

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
%springDat.nep = nep;

