load('MPN.mat')
load([MPN 'obI.mat']);

%% Get skeletons
skelList = 125;
minTip = 2;
minSpace = 2;

    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,skelList);
    load(fileName)
    %%
    skel = cellStruct.arbor;
    node2subs = skel.nodes.node2subs;
    bones = skel.branches;
    
    edges = cat(1,bones.edges, skel.bridges);
    nodes = 1:size(node2subs,1); % will go wrong if gaps
    
    nodePos = node2subs(nodes,:);
    
    voxelScale = obI.em.dsRes * 2; %adjusted for skeleton down samp
    for d = 1:size(nodePos,2)
        nodePos(:,d) = nodePos(:,d) * voxelScale(d);
    end
    
    nep.nodePos = nodePos;
    nep.nodes = nodes;
    nep.edges = edges;
    nep.cellID = skelList(1);
    


    %% condition skeleton
        
   % showNep(nep)
    clf
    subplot(1,3,1)
    nep = uniteSkel(nep);%% fix breaks in skeleton
    nep = groupNepEdges(nep); %% group edges into bones
    nep = edgeGroup2bones(nep);
    nep = bonesVsSpurs(nep,minTip); %% sort bones and spurs
    nep = simpleNep(nep,minSpace);
    showNepBones(nep,2);
    pause(.01)

        subplot(1,3,2)

    %%redo without spurs
    nep.edges = cat(1,nep.bones.edges);
    nep = groupNepEdges(nep);
    nep = edgeGroup2bones(nep);
    showNepBones(nep,2);
    pause(.01)
    
        subplot(1,3,3)

    %%Simplify skeleton
    nep = cleanBones(nep);
    nep = uniteSkel(nep);%% fix breaks in skeleton
    showNepBones(nep);
    pause(.01)
    
    if 0
        nepFileName = sprintf('%snep\\skelNep%d.mat',MPN,skelList);
        save( nepFileName,'nep')
    end
    
    
    
    
    
    
    