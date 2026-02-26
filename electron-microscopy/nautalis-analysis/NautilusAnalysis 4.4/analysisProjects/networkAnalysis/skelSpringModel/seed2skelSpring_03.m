



%% Get data
loadData = 1;
if loadData
    clear all
    %MPN = GetMyDir;
    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    cellList = unique([plusOne; useList.preList(:); useList.postList(:)]);
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end
binaryMat = 0;


%load([springRes 'turkey_plusOne125.mat'])


%% Get skeletons
skelList = 125;
minTip = 10;
minSpace = 10;

for i = 1 : length(skelList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,skelList(i));
    load(fileName)
    %%
    skel = cellStruct.arbor;
    node2subs = skel.nodes.node2subs;
    bones = skel.branches;
    
    edges = cat(1,bones.edges, skel.bridges);
    nodes = unique(edges(:)); % will go wrong if gaps
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

%% Merge skel, cell and syn

nepCSS = mergeSkelCellSyn(skels,nepCell,nepSyn);
showNepCSS(nepCSS);

nepName = 'nepTest_128.mat';
nepDir = 'D:\LGNs1\Analysis\LIN\nepDir\';
%%%% save([nepDir nepName],'nepCSS');


%% Filter nodes
%nodes = useList.nodes;
seedFilt = 108;
seedPos = find(nodes==seedFilt);
con = useList.con;
symCon = con + con';
[cellGroup] = segmentCon(symCon);
seedGroup = cellGroup(seedPos);
nodes = nodes(cellGroup == seedGroup);

nodeType = nodes*0;
for i = 1:length(nodes)
    nodeType(i) = useList.nodeType(useList.nodes==nodes(i));
end
nodeIDs = nodes;




nodeNum = length(nodeIDs);
lookUpID(nodeIDs+1) = 1:length(nodeIDs);


allEdges(:,3) = 1;
allEdges(allEdges(:,1) == plusOne,3) = 1;
allEdges(allEdges(:,2) == plusOne,3) = 1;

%% Set color
nodeCol = zeros(nodeNum,3);

if 1
    
    [nodeCol, useCol] = getList_seedColor(seedList,nodeIDs) ;
    nodeCol = nodeCol * .8;
    nodeCol(sum(nodeCol,2) == 0,:) = .3;
    
    seedCol = [1 0 0; 0 1 0; .5 0 1 ; 0 .5 1];
    for s = 1:length(seedList)
        targ = find(nodeIDs == seedList(s));
        nodeType(targ) = nodeType(targ)+2;
        nodeCol(targ,:) = seedCol(s,:);
    end
    
    showCol = permute(useCol.col,[1 3 2])
    
    image(uint8(showCol*256))
    useCol.comb{10}
    
    targ = find(nodeIDs == plusOne);
    nodeCol(targ,:) = [10 10 10]
    
end


%%Set color according to attribute
if 0
    %[colorList cellCol] = getAttributes(obI);
    
    load([MPN 'cb2d.mat'])
    colorList = cb2d.IDs;
    nodeCol = nodeCol + .2;
    
    colPropRaw = cb2d.areaUM;
    
    nodeProp = [];nodePropRef = [];
    for i = 1:length(nodeIDs)
        targ = find(colorList == nodeIDs(i));
        if ~isempty(targ)
            if length(targ)>1
                'too many targets'
                colorList(targ)
            end
            nodeProp= [nodeProp colPropRaw(targ(1))];
            nodePropRef = [nodePropRef i];
            
        end
    end
    
    colProp = nodeProp-min(nodeProp)+1;
    colProp = ceil(colProp/max(colProp) * 100);
    
    colTable = jet(100);
    cellCol = colTable(colProp,:);
    
    nodeCol(nodePropRef,:) = cellCol;
    nodeType(nodePropRef) = nodeType(nodePropRef) + 2;
    
    cat(2,nodeIDs', nodeCol);
end

%%Label specific cells
if 0
    crossoverAxons = [2032	2033	2034	2035]
    gotList = getList_giantBoutons(MPN);
    labelCells = gotList;
    for i = 1:length(labelCells)
        nodeCol(find(nodeIDs==labelCells(i)),:) = nodeCol(find(nodeIDs==labelCells(i)),:) + .5;
    end
end


nodeCol(nodeCol>1) = 1;
nodeCol(nodeCol<0) = 0;

%% Create springDat (need nodeCol, nodeIDs,

springIn.nodeIDs = nodeIDs;
springIn.allEdges = allEdges;
springIn.allWeights = allEdges ;
springIn.nodeCol = nodeCol;
springIn.nodeType = nodeType;
springIn.seedList = seedList;

springDat = skelSpringParameters_01(springIn);

if binaryMat
    springDat.edges.ew = springDat.edges.ew>0;
end

%% Run springs
movDir  = [springDir 'seedMovie_?\'];
if ~exist(movDir,'dir'), mkdir(movDir), end

for rerun = 1: 1
    allResults{rerun} = runSkelSprings(springDat);
end

%%  print eps
if 0
    springDir = 'D:\LGNs1\Analysis\springDat\skelSpring\';
    if ~exist(springDir,'dir'), mkdir(springDir), end
    
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    tag = 'net125_white';
    
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
end



%% Save result
if 0
    
    springRes = 'D:\LGNs1\Analysis\springDat\results\';
    if ~exist(springRes,'dir'), mkdir(springRes), end
    
    result = allResults{end};
    save([springRes 'turkey_plusOne125.mat'],'result')
    
    
end






