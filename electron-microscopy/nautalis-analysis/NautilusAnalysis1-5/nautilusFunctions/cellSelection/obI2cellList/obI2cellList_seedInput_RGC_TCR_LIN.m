function[useList] = obI2cellList_seedInput_RGC_TCR_LIN(obI,seedList);

%%Select a list of pre and post synaptic cells based on obI (object
%%information), a seed input, and a scripted set of filters.


%load([MPN 'obI.mat'])

%% variables
%seedList = [108]
crossoverAxons = [2032	2033	2034	2035];

%% Get cells
tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
linList = obI.nameProps.cellNum(obI.nameProps.lin);
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

preList = [];
for i = 1:length(seedList)
    isPost = edges(:,1) == seedList(i);
    preList{i} = unique(edges(isPost,2));
end
%axList =    setxor(preList{1},preList{2});  %dont use a
axList =    cat(1,preList{:});  %dont use a
% 
axList = intersect(unique(axList),[rgcList linList]);
% axList = axList((axList>=1000) & (axList<10000));
axList = axList(axList>0);
% axList = setdiff(axList,noSkel);
% disp(sprintf('Results calculated without %d',noSkel));

% axList = setdiff(axList, crossoverAxons)
% axList = axList((axList>=1000) & (axList<5000));
% 
% axList = axList(axList<2000);
% axList = [axList crossoverAxons]

postList = [];
for i = 1:length(axList)
    isPre = edges(:,2) == axList(i);
    postList = [postList; edges(isPre,1)];
end
%cellList = unique(postList);
cellList = intersect(unique(postList),[tcrList linList]);
cellList = cellList(cellList>0);
cellList = cellList((cellList>0) & (cellList<1000));
cellList = [seedList setdiff(cellList,seedList)];

%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end


%% Create outpute structure

useList.seedList = seedList;
useList.preList = axList;
useList.postList = cellList;
useList.con = con;




