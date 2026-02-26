

load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
    aspect =  [1.2267 1.0667 1]

    tracedCells = [];
    showTraced = intersect(tracedCells, [obI.cell.name])
    
%% make OBJ
cellList =  num2cell(showTraced(showTraced>9),[1 length(showTraced)]);
renderOb = 0;
%cellList = [10	129	162	170]
objDir = [MPN '\objFiles\'];

objDir = [MPN 'obFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end

cellSubs = names2Subs(obI,dsObj,cellList);
downSamp = 4;
cellDat.cellList = cellList;
for i = 1:length(cellSubs)
    sub = cellSubs{i};
    for a = 1:3
        sub(:,a) = sub(:,a)* aspect(a);
    end
    
    
    obName = cellList{i};
    if ~ischar(obName),obName = num2str(obName);end
    if ~isempty(sub)
    smallSub = shrinkSub(sub,downSamp);
    tic
    fv = subVolFV(smallSub,[],renderOb);
    fileName = sprintf('%sdSamp%d_%s.obj',objDir,downSamp,obName);
    vertface2obj(fv.vertices,fv.faces,fileName,obName);
    toc
    cellDat(i).subs = sub;
    cellDat(i).fv = fv;
    end
    
end

