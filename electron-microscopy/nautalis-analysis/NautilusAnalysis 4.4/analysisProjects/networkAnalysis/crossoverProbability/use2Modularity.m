% calculate the modularity of axons where their group is defined by
% connectivity to a seed cell and edges in the adjacency matrix are defined
% by axons converging on the same target cell.  Self edges are excluded
% useList = obI2cellList_seedInput_RGC_TCR(obI,seedCells)

function[Q]  = use2Modularity(useList);



% MPN = GetMyDir
% load([MPN 'obI.mat'])
% seedCells = [108 201];
% useList = obI2cellList_seedInput_RGC_TCR(obI,seedCells);
conPref = seedPreferences(useList.seedList,useList);

rawCon = useList.con;
[a isSeed] = intersect(useList.postList,useList.seedList);
notSeed = setdiff(1:length(useList.postList),isSeed);
seedCon = rawCon(:,isSeed);
con = rawCon(:,notSeed);


%% get links
links = zeros(size(con,1),size(con,1));
ingroup = links;
for y = 1:size(con,1)
    for x = 1:size(con,1)
        links(y,x) = sum((con(y,:)>0) & (con(x,:)>0));
        if y == x
           links(y,x) = 0;
        end
        ingroup(y,x) = (seedCon(y,1)>0) == (seedCon(x,1)>0);
    end
end

% 
 image(links*20)
% image(ingroup*30)

%%

adj = links;
modules = {[find(seedCon(:,1)>0)] [find(seedCon(:,2)>0)]};
Q=modularity_metric(modules,adj);






