function[sm] = makeSkeletonForSM(sm);

%% Make skeleton
%%To Do: Improve seed, Improve conMat Bridging. Improve arbor format,
%%electrotonic with rads. update simplify
load('MPN.mat');
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

rawObjectSubs = getCellSubs(obI,dsObj,sm.cid);
rawSeed   = double(sm.cell.anchor);
fixSeed = rawSeed .* obI.em.res ./obI.em.dsRes/1000;
dsamp = [1 1 1];
if (size(rawObjectSubs,1)>0)
    
    objectSubs = double(downSampSub(rawObjectSubs,dsamp));
    
    if 1
        oneEnd = rawObjectSubs(rawObjectSubs(:,3) == max(rawObjectSubs(:,3)),:);
        rawSeed = median(rawSeed,1);
        seedSub = rawSeed./dsamp;
    end
    if isempty(rawSeed)
        seedSub = median(objectSubs,1);
    else
        seedSub = round(fixSeed./dsamp);
    end
    
    
    cellStruct = subs2arbor(objectSubs,seedSub);
    cellStruct.voxSize = obI.em.dsRes.*dsamp;
    cellStruct.arbor.voxSize = obI.em.dsRes.*dsamp;

    
    stopTime = clock;
    
else
    'no subs found'
end

arbor = cellStruct.arbor
%showRadSurf(arbor.nodes.pos,arbor.edges, cellStruct.arbor.nodes.rad)


%% reprocess skeleton
minTip = 0;
minSpace = 0;
%%


arbor = cellStruct.arbor;
edges = arbor.edges;
nodes = arbor.nodes.nodes; % will go wrong if gaps

nodePos = arbor.nodes.pos;
clear nep
nep.pos = nodePos * arbor.voxSize(1);
nep.nodes = arbor.nodes.nodes;
nep.edges = edges;
nep.edgeRad  = arbor.edgeProps.rad * arbor.voxSize(1);
nep.nodeRad = arbor.nodes.rad * arbor.voxSize(1);

%showRadSurf(nep.pos,nep.edges, nep.nodeRad)


%% condition skeleton
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
nep = groupNepEdges(nep);
%showRadSurf(nep.pos,nep.edges, nep.nodeRad)

nep = edgeGroup2bones(nep)
%showRadSurf(nep.pos,nep.edges, nep.nodeRad)

sm.nep = nep;
sm.arbor = arbor;



