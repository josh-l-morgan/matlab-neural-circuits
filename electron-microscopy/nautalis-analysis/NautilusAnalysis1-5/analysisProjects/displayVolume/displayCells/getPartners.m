function[cellList] = getPartners(obI,sourceCells,prePostAll)


if ~exist('prePostAll','var'),prePostAll = 'all';end
targCells = unique(sourceCells(sourceCells>0));

postCellList = [];
preCellList = [];

for i = 1: length(targCells)
    postHits =  find(obI.nameProps.edges(:,2) == targCells(i));
    postCellList = cat(1,postCellList,obI.nameProps.edges(postHits,1)); 
    preHits =  find(obI.nameProps.edges(:,1) == targCells(i));
    preCellList = cat(1,preCellList,obI.nameProps.edges(preHits,2)); 
        
end


if strcmp(lower(prePostAll),'pre')
    cellList = preCellList;
elseif strcmp(lower(prePostAll),'post')
    cellList = postCellList;
else
    cellList = cat(1,preCellList,postCellList);
end


cellList = unique(cellList(cellList>0));
