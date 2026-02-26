



%% Get data
loadData = 1;
if loadData
    clear all
    %MPN = GetMyDir;
    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903 125];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
   % useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
   useList = obI2nodes_rtl_plus1(obI,plusOne);

    cellList = unique([plusOne; useList.preList(:); useList.postList(:);forceNodes]);
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end
binaryMat = 0;


%load([springRes 'turkey_plusOne125.mat'])


%% Get skeletons
skelList = 125;
minTip = 2;
minSpace = 2;

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
    voxelScale = obI.em.dsRes * 2; %adjusted for skeleton down samp
    nodePos = scaleSubs(nodePos,voxelScale);
    nep.nodePos = nodePos;    
    
    nep.nodes = nodes;
    nep.edges = edges;
    


    %% condition skeleton
        
    showNep(nep)

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
WPN = 'D:\LGNs1\Analysis\'

nepName = 'nepTest_125c.mat';
nepDir = [WPN 'LIN\nepDir\'];
 save([nepDir nepName],'nepCSS');




