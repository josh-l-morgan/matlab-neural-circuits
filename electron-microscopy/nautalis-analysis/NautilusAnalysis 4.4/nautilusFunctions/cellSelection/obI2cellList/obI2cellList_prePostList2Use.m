function[useList] = obI2cellList_prePostList2Use(obI,axList,cellList);

%%Make list of all traced thalamocortical cells and the rgcs that innervate
%%them
if ~exist('seedList','var')
    seedList = [];
end


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

useList.seedList = seedList(:)';
useList.preList = axList(:)';
useList.postList = cellList(:)';
useList.con = con;




