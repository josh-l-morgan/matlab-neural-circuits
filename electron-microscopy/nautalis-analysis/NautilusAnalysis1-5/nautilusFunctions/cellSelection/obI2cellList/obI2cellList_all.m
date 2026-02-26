function[useList] = obI2cellList_all(obI);

%%Select a list of pre and post synaptic cells based on obI (object
%%information), a seed input, and a scripted set of filters.

if ~exist('seedList','var')
    seedList = [];
end
axList =  obI.cell.name(obI.cell.name>0);
cellList =  obI.cell.name(obI.cell.name>0);
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

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




