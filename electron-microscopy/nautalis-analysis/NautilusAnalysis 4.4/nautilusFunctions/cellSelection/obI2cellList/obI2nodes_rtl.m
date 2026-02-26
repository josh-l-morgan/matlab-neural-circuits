function[useList] = obI2nodes_rtl(obI);


%% Get cells



tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
linList = obI.nameProps.cellNum(obI.nameProps.lin);
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);


%% get nodes
nodes = unique(obI.cell.name);
nodeType = nodes * 0;
for i = 1:length(nodes)
    targ = find(obI.cell.name == nodes(i),1);
    if obI.nameProps.rgc(obI.cell.mainObID(targ))
        nodeType(i) = 1;
    elseif obI.nameProps.tcr(obI.cell.mainObID(targ))
        nodeType(i) = 2;
    elseif obI.nameProps.lin(obI.cell.mainObID(targ))
        nodeType(i) = 3;
    else
        nodeType(i) = 4;
    end
    
end

useNodes = (nodes>0) & (nodeType>0) & (nodeType < 30);
nodes = nodes(useNodes);
nodeType = nodeType(useNodes);

nodes = [nodes];
nodeType = [nodeType];

%% graph

con = zeros(length(nodes),length(nodes));
for i = 1:length(nodes)
    for p = 1:length(nodes)
        con(i,p) =  con(i,p)+ sum( (edges(:,1) == nodes(p)) & (edges(:,2) == nodes(i)));
        con(i,p) =  con(i,p)+ sum( (edges(:,1) == nodes(i)) & (edges(:,2) == nodes(p)));

    end
end

%% Create outpute structure

useList.seedList = nodes;
useList.preList = nodes;
useList.postList = nodes;
useList.nodes = nodes;
useList.nodeType = nodeType;
useList.con = con;


% 
% preList = [];
% for i = 1:length(seedList)
%     isPost = edges(:,1) == seedList(i);
%     preList{i} = unique(edges(isPost,2));
% end
% axList = cat(1,preList{:});  %dont use a
% axList = intersect(unique(axList),[rgcList linList]);
% axList = axList(axList>0);
% 
% postList = [];
% for i = 1:length(axList)
%     isPre = edges(:,2) == axList(i);
%     postList = [postList; edges(isPre,1)];
% end
% cellList = intersect(unique(postList),[tcrList linList]);
% cellList = cellList(cellList>0);
% cellList = [seedList setdiff(cellList,seedList)];










