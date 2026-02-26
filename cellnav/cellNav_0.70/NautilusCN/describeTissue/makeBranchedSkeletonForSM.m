function[sm] = makeBranchedSkeletonForSM(sm);

%% Make skeleton
%%To Do: Improve seed, Improve conMat Bridging. Improve arbor format,
%%electrotonic with rads. update simplify
load('MPN.mat');
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

dsamp = [1 1 1];
cellSub = obj2conBranch(obI,dsObj,sm.cid);

if ~isempty(cellSub)

%%Soma
cellSub.soma = obIobSoma(obI,dsObj,sm.cid);
renderCon(cellSub.subs + 1,[],[0 0 1])
cellSub.isSoma = zeros(size(cellSub.subs,1),1);
if ~isempty(cellSub.soma.subs)
    if 1;%sum(cellSub.soma.subs>0) == 3;
    renderCon(cellSub.soma.subs,[],[1 0 0])
    
     
    somY = zeros(max(cellSub.subs(:,1)),1);
    somX = zeros(max(cellSub.subs(:,2)),1);
    somZ = zeros(max(cellSub.subs(:,3)),1);
    somY(cellSub.soma.subs(:,1)) = 1;
    somX(cellSub.soma.subs(:,2)) = 1;
    somZ(cellSub.soma.subs(:,3)) = 1;
    isY = somY(cellSub.subs(:,1));
    isX = somX(cellSub.subs(:,2));
    isZ = somZ(cellSub.subs(:,3));
    cellSub.isSoma = (isY + isX + isZ) == 3;
    clf
    renderCon(cellSub.subs(cellSub.isSoma==1,:),[],[0 0 1])
    renderCon(cellSub.subs(~cellSub.isSoma==1,:),[],[1 0 0])
    end
end




%%Seed
anc2sub = (obI.em.res / 1000)./ obI.em.dsRes ;
obI.em.dsRes.*dsamp;
targ = find(obI.cell.name == sm.cid);
if obI.cell.seedID{targ}
    anch = double(obI.cell.anchors(targ,[2 1 3]));
    pos = anch .* anc2sub;
else
    pos = cellSub.soma.center;
end

dists = getDist(cellSub.subs,pos);
seedInd = find(dists==min(dists),1);
seedSub = cellSub.subs(seedInd,:);
cellSub.seedSub = seedSub;
cellSub.seedInd = seedInd;

cellSub = parseCellSubs(cellSub);
%renderCon(cellSub.subs,cellSub.conMat)
cellSub = bridgeConMat(cellSub);


%objectSubs = double(downSampSub(rawObjectSubs,dsamp));
cellStruct = subs2arbor(cellSub);
cellStruct.voxSize = obI.em.dsRes.*dsamp;
cellStruct.arbor.voxSize = obI.em.dsRes.*dsamp;
arbor = cellStruct.arbor;

%%Show Result
clf
renderCon(arbor.vox.subs,[],[1 0 0],.1)
showRadSurf(arbor.nodes.pos,arbor.edges, cellStruct.arbor.nodes.rad,[0 0 1],.1)
showRadSurf(arbor.nodes.pos,arbor.edges, cellStruct.arbor.nodes.rad*0+.1,[0 1 0],1)

%% reprocess skeleton
arbor = cellStruct.arbor;
seedPos = arbor.vox.subs(arbor.skel.surfSeed,:) ;


edges = arbor.edges;
nodes = arbor.nodes.nodes; % will go wrong if gaps

nodePos = arbor.nodes.pos;
dists = getDist(nodePos,seedPos);

clear nep
nep.seedNode = find(dists==min(dists));
nep.pos = nodePos * arbor.voxSize(1);
nep.nodes = arbor.nodes.nodes;
nep.edges = edges;
nep.edgeRad  = arbor.edgeProps.rad * arbor.voxSize(1);
nep.nodeRad = arbor.nodes.rad * arbor.voxSize(1);

nep = groupNepEdges(nep);
%showRadSurf(nep.pos,nep.edges, nep.nodeRad)

nep = edgeGroup2bones(nep)
nep = smoothRad(nep,.2);
nep.fv = arbor.vox.fv;
nep.fv.vertices = nep.fv.vertices * arbor.voxSize(1);

clf
renderFV(nep.fv,[1 0 0],.2)
%showRadSurf(nep.pos,nep.edges, nep.nodeRad,[0 1 0],.4)
showRadSurf(nep.pos,nep.edges, nep.meanNodeRad, [ 0 0 1],.2)
showRadSurf(nep.pos,nep.edges, nep.meanNodeRad * 0 + .01, [ 0 1 0],1)


sm.nep = nep;
sm.arbor = arbor;

%% unused analysis
%showRadSurf(nep.pos,nep.edges, nep.nodeRad)
%%condition skeleton
%{
showNep(nep)

%nep = uniteSkel(nep);%% fix breaks in skeleton
nep = groupNepEdges(nep); %% group edges into bones
nep = edgeGroup2bones(nep);
%nep = bonesVsSpurs(nep,minTip); %% sort bones and spurs

%%redo without spurs
% nep.edges = cat(1,nep.bones.edges);
% nep = groupNepEdges(nep);
% nep = edgeGroup2bones(nep);
% showNepBones(nep,2);

%%Simplify skeleton
nep = simpleNep(nep,minSpace);
showNepBones(nep)
nep = cleanBones(nep);
showNepBones(nep);
%}

else
    sm.nep = [];
    sm.arbor = [];
    
end



