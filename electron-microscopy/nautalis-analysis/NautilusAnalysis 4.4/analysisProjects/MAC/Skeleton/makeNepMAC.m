


%% Get data
loadData = 1;
if loadData
    clear all
    %MPN = GetMyDir;
    MPN = 'Z:\joshm\LGNs1\Exports\MAC_export\MAC_merge_mat\';
    WPN = MPN;

    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903 125];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
   % useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);

    allEdges = obI.nameProps.edges(:,[2 1]);
    
end
binaryMat = 0;

postList = allEdges(:,2);
postList = unique(postList(postList>0));
%load([springRes 'turkey_plusOne125.mat'])

cells = obI.cell.name(obI.cell.isCell>0);
cellList = intersect(postList,cells);

%% Get skeletons
skelList = cellList;
minTip = 0;
minSpace = 1;

for i = 1 : length(skelList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,skelList(i));
    load(fileName)
    %%
   skel = cellStruct.arbor;
    node2subs = skel.nodes.node2subs;
    bones = skel.branches;
    
    edges = cat(1,bones.edges, skel.bridges);
    nodes = 1:size(node2subs,1); % will go wrong if gaps
    
    nodePos = node2subs(nodes,:);
    voxelScale = obI.em.dsRes; %adjusted for skeleton down samp
    nodePos = scaleNodes(nodePos,voxelScale);
    nep.nodePos = nodePos;    
    
    nep.nodes = nodes;
    nep.edges = edges;
    


    %% condition skeleton
        
    showNep(nep)

   % nep = uniteSkel(nep);%% fix breaks in skeleton
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
    
    skels(i).id = skelList(i);
    skels(i).nep = nep;
    
end

%% Load synapses

nepSyn = obISyn2nep(obI);
nepSyn = linkSyn2Skel(nepSyn,skels);

%% get cell starting positions
[cellID cellProp] = getList_cellPositions;
nepCell.nodes = [cellList];
for i = 1:length(cellList);
    foundPos = cellProp(cellID == cellList(i),:);
    if ~isempty(foundPos)
        cellPos(i,:) = foundPos;
    else
%         foundPos = find(obI.cell.name == cellList(i));
%         
%         cellPos(i,:) = obI.cell(foundPos).anchors
        cellPos(i,:) = mean(cellProp,1);
        'missing position'
        
    end
    
end
nepCell.nodePos = cellPos;
nepCell.nodeParent = cellList;

%% Merge skel, cell and syn

nepCSS = mergeSkelCellSyn(skels,nepCell,nepSyn);
nepCSS.obI = obI;
showNepCSS(nepCSS);

nepName = 'nepMAC.mat';
nepDir = [WPN '\nepDir\'];
if ~exist(nepDir,'dir'),mkdir(nepDir); end
 save([nepDir nepName],'nepCSS');




