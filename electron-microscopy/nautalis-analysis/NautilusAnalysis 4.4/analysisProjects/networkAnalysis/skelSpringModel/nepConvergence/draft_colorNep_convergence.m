function[nep] = colorNep(nep)



%% Test color

%%nodes
%%%%%%% cell syn skel %%%%%%%
nodeNum = length(nep.nodes);
nep.nodeCol = ones(nodeNum,3)*.2;
nep.nodePushes = ones(nodeNum,1);
nep.nodeIsPushed = ones(nodeNum,1);
nep.nodeMass = ones(nodeNum,1);
nep.nodeMove = ones(nodeNum,1);

nodeTypeList = unique(nep.nodeType);
colList = [.3 .3 .3 ;0 1 1; 0 .6 .6;  .8 .5 .5; .5 .5 .5 ; .5 .5 .5];
pushes = [1 .2 .1 0 0];
isPushed = [1 0 0 0.01 0  0];
mass = [1 1 1 1 1 1];
usePos = [1 1 1 1 1]>0;
nodeMove = [ 1 0 0 0 0 ]>0;
nodeDisplay = [0 0 0 1 1];
markerList = {'o','^','s','d','p','h','.','<','>','+','*','x'}
nodeMarker = {'^', 'o', '.'}
for i = 1:length(nodeTypeList)
    isType = find(nep.nodeType == nodeTypeList(i));
    fullColor = repmat(colList(nodeTypeList(i),:),[length(isType) 1]);
    nep.nodeCol(isType,:) = fullColor;
    nep.nodePushes(isType) = pushes(nodeTypeList(i));
    nep.nodeIsPushed(isType) = isPushed(nodeTypeList(i));
    nep.nodeMass(isType) = mass(nodeTypeList(i));
    nep.usePos(isType) = usePos(nodeTypeList(i));
    nep.nodeMove(isType) = nodeMove(nodeTypeList(i));
    nep.nodeDisplay(isType) = nodeDisplay(nodeTypeList(i));
    nep.nodeMarker(isType) = nodeMarker(nodeTypeList(i));

end


if 1
    
    conTarg = 125;
    edges = nep.obI.nameProps.edges(:,1:2);
    nodeNames = nep.nodeName;
    
    for i = 1:length(nodeNames)
       postEdges = edges((edges(:,1)==nodeNames(i)) & (edges(:,2)==conTarg),:)
       preEdges = edges((edges(:,2)==nodeNames(i)) & (edges(:,1)==conTarg),:)
       preCon(i) = size(preEdges,1);
       postCon(i) = size(postEdges,1);
       if ~nodeNames(i)
           preCon(i) = 0;
           postCon(i) = 0;
       end
    end
    
    
    
    [nodeCol useCol] = getList_seedColor(1,nep.nodeName) ;
    nodeCol = nodeCol * .8;
    useCol = sum(nodeCol,2)>0;
    nep.nodeCol(useCol,:) = nodeCol(useCol,:);
    nep.nodeDisplay(useCol) = 1;
    foundNodes = find(useCol);
    
    
end


if 0 %use getList_seeedColor
    [nodeCol useCol] = getList_seedColor(1,nep.nodeName) ;
    nodeCol = nodeCol * .8;
    useCol = sum(nodeCol,2)>0;
    nep.nodeCol(useCol,:) = nodeCol(useCol,:);
    nep.nodeDisplay(useCol) = 1;
    foundNodes = find(useCol);
end

if 1 %Assign cell type markers
    
    [checkIDs checkProp]  = getList_cellTypes;
    findCells = find(nep.nodeType==1);
    cellList = nep.nodeParent(findCells);
    [int a b] = intersect(checkIDs, cellList);
    cellType = checkProp(a);
    typeMarkers = {'.','^','o','s','p','h'};
    cellMarkers = typeMarkers(cellType+1);
    nep.nodeMarker(findCells) = cellMarkers;
end

%% edges
edgeNum = size(nep.edges,1);
nep.edgeCol = ones(edgeNum,3)*.1;
nep.edgeWeights = ones(edgeNum,1);

%%Default
% pre post skel cell2skel
edgeTypeList = unique(nep.edgeType);
nep.edgeWidth = ones(edgeNum,1);
colList = [.2 .2 .2; .2 .2 .2; 0 .5 .5 ; .2 1 1; .5 .5 .5 ; .5 .5 .5];
widthList = [ 1 1 2 1 1 1];
weightList = [10^10 1 1 1];
edgeDisplay = [0 0 1 1 1];
for i = 1:length(edgeTypeList)
    isType = find(nep.edgeType == edgeTypeList(i));
    fullColor = repmat(colList(edgeTypeList(i),:),[length(isType) 1]);
    nep.edgeCol(isType,:) = fullColor;
    nep.edgeWidth(isType,:) = widthList(edgeTypeList(i));
    nep.edgeWeights(isType,:) = weightList(edgeTypeList(i));
    nep.edgeDisplay(isType,:) = edgeDisplay(edgeTypeList(i));
end


if 1 % color edges according to pre
    
    isPre = nep.edgeType == 1; %find pre edges
    nep.edgeCol(isPre,:) = nep.nodeCol(nep.edges(isPre,1),:); %set pre edges to  color of node of pre
    
    isPost = nep.edgeType == 2;
    nep.edgeCol(isPost,:) = nep.nodeCol(nep.edges(isPost,2),:);
    nep.nodeCol(nep.edges(isPost,1),:) = nep.nodeCol(nep.edges(isPost,2),:); % color synapse according to pre
   
nep.edgeDisplay(isPost) = nep.edgeDisplay(nep.edges(isPost,2));
   nep.nodeDisplay(nep.edges(isPost,1)) = nep.nodeDisplay(nep.edges(isPost,2)); % display synapse according to pre
   
   
%     nep.edgeCol = nep.nodeCol(nep.edges(:,1),:);
%     nep.edgeWidth = ones(edgeNum,1);
end


if 1 % Display edges connected to displayed nodes
   
    %foundNodes = find((nep.nodeDisplay == 1) & (nep.nodeType == 1));
%     [i a b] = intersect(nep.edges(:,1),foundNodes);
   % nep.edgeDisplay(a) = 1;
    [i a b] = intersect(nep.edges(:,2), foundNodes);
    nep.edgeDisplay(a) = 1;
        
end


if 0 % color pre post edges of 125
    from125 = find(nep.nodeParent == 125);
    nep.edgeCol(find(ismember(nep.edges(:,1),from125) & (nep.edgeType == 1)),1) = 1;
    nep.edgeCol(find(ismember(nep.edges(:,2),from125) & (nep.edgeType == 2)),2) = 1;
end

%% Color Target 

if 0 % color pre post nodes of 125
    targCell = 125;
    fromTarg = find(nep.nodeParent == targCell); %list of nodes belonging to target cell
    preTarg = find(ismember(nep.edges(:,1),fromTarg) & (nep.edgeType == 1)); %targ is pre
    postTarg = find(ismember(nep.edges(:,2),fromTarg) & (nep.edgeType == 2)); %edge where targ is post
    
    %nep.edgeCol(preTarg,:) = 1;
    %nep.edgeCol(postTarg,:) = 1;

    preNodes = nep.edges(preTarg,2);
    postNodes = nep.edges(postTarg,1);
   
    nep.nodeCol(preNodes,:) = repmat([0 1 1],[length(preNodes) 1]);
    nep.nodeCol(postNodes,:) = repmat([1 1 0],[length(postNodes) 1]);;
   %nep.nodeCol(fromTarg,:) = 1;
end




%% Change type contrasts

nep.nodeCol(nep.nodeType == 1,:) = nep.nodeCol(nep.nodeType == 1,:) * .3;
nep.nodeCol(nep.nodeType == 2,:) = nep.nodeCol(nep.nodeType == 2,:) * 1;
nep.nodeCol(nep.nodeType == 3,:) = nep.nodeCol(nep.nodeType == 3,:) * 1;
nep.edgeCol(nep.edgeType == 1,:) = nep.edgeCol(nep.edgeType == 1,:) * .3;
nep.edgeCol(nep.edgeType == 2,:) = nep.edgeCol(nep.edgeType == 2,:) * .3;
nep.edgeCol(nep.edgeType == 3,:) = nep.edgeCol(nep.edgeType == 3,:) * 1;
nep.edgeCol(nep.edgeType == 4,:) = nep.edgeCol(nep.edgeType == 4,:) * 1;
nep.edgeCol(nep.edgeType == 5,:) = nep.edgeCol(nep.edgeType == 5,:) * 1;



%{


if 0
    
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

nep.nodeCol = nodeCol;

%}