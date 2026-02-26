%{
MPN = GetMyDir;
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
%}
voxelScale = [1 1 1]; %adjusted for skeleton down samp

%% Get subs

objectSubs = double(dsObj(1).subs);
scatter(objectSubs(:,1),objectSubs(:,2),'.')
seedSub = objectSubs(1,:);

%% Skeletonize
cellStruct = subs2arbor(objectSubs,seedSub);

%% convert to NEP
skel = cellStruct.arbor;
node2subs = skel.nodes.node2subs;
bones = skel.branches;

edges = cat(1,bones.edges, skel.bridges);
nodes = unique(edges(:)); % will go wrong if gaps
nodePos = node2subs(nodes,:);
nodePos = scaleSubs(nodePos,voxelScale);

nep.nodePos = nodePos;
nep.nodes = nodes;
nep.edges = edges;



%% condition skeleton

showNep(nep)
minTip = 50;
minSpace = 5;
nep = uniteSkel(nep);%% fix breaks in skeleton
nep = groupNepEdges(nep); %% group edges into bones
nep = edgeGroup2bones(nep);
nep = bonesVsSpurs(nep,minTip); %% sort bones and spurs

%%redo without spurs
nep.edges = cat(1,nep.bones.edges);
nep = groupNepEdges(nep);
nep = edgeGroup2bones(nep);
showNepBones(nep,2);

%%Simplify skeleton
nep = simpleNep(nep,minSpace);
showNepBones(nep)
nep = cleanBones(nep);
showNepBones(nep);
