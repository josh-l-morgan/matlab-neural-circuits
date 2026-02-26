function[useLists] = obI2cellList(obI);


%load([MPN 'obI.mat'])

%% variables
seedList = [108]
crossoverAxons = [2032	2033	2034	2035]
noSkel = [2014 1026]

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

axList = intersect(unique(axList),rgcList);
axList = axList((axList>=1000) & (axList<5000));
axList = setdiff(axList,noSkel);
disp(sprintf('Results calculated without %d',noSkel));

% axList = setdiff(axList, crossoverAxons)
% axList = axList((axList>=1000) & (axList<5000));
% 
% axList = axList(axList<2000);
% axList = [axList crossoverAxons]

postList = [];
for i = 1:length(axList)
    isPre = edges(:,2) == axList(i);
    postList = [postList; edges(isPre,:)];
end
cellList = intersect(unique(postList),tcrList);
cellList = cellList((cellList>0) & (cellList<1000));

%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end
