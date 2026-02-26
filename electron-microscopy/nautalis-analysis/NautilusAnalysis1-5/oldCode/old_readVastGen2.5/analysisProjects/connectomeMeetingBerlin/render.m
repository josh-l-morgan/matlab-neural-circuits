



%% make OBJ
cellList = [108 201];
cellList = num2cell(cellList);
renderOb = 1;
%cellList = [10	129	162	170]
%objDir = ['D:\LGNs1\Segmentation\VAST\S4\cellLibrary\objFiles\'];

%objDir = [MPN 'obFiles\']
%if ~exist(objDir,'dir'),mkdir(objDir),end

cellSubs = names2Subs(obI,dsObj,cellList);
downSamp = 16;
cellDat.cellList = cellList;
for i = 1:length(cellSubs)
    sub = cellSubs{i};
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


%%

for i = 1: length(cellDat)
    fv = cellDat(i).fv;
p = patch(fv);
view(30,-15);
axis vis3d;
colormap copper
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
end


end