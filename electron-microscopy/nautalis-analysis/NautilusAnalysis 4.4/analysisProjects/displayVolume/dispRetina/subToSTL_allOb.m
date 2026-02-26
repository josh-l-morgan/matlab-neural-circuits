clear all
load('MPN.mat')
load('WPN.mat')
%dSamp = [

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

targCell = 6;
allEdges = obI.nameProps.edges;
cellNames  = unique(obI.nameProps.cellNum);
cellList = [1:length(dsObj)];
% 
% colMap = hsv(256);
% colNum = length(cellList)-1;
% rainBow = colMap(floor([1:colNum] * 255/(colNum)),:);
% col = [rainBow; [1 1 1]]
%col = [1 1 1;repmat([1 0 0],[length(preTarg) 1]);repmat([0 1 0],[length(postTarg) 1])]
col = hsv(length(cellList));

renderOb = 0;
tag = 'testCrop';
objDir = [WPN 'stlFiles\wafovbas\']

if ~exist(objDir,'dir'),mkdir(objDir),end

downSamp = 1;

target = [908 490 555.5]*2; % ( Y X Z)
rad = 300;
%crop = [target- rad; target + rad];
%crop     = [2000 400 600;          2500 600 700];
clf
for i = 1:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
    sub = dsObj(i).subs;
    obName = cellList(i);
    if iscell(obName); obName = obName{1};end
    if exist('crop','var')
        useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
            (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
            (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
        sub = sub(useSub,:);
        
    end
    smallSub = shrinkSub(sub,downSamp);
    tic
    if isempty(smallSub)
%        disp(sprintf('no points on %d',cellList{i}))
    else
    fv = subVolFV(smallSub,[],renderOb);
    renderFV(fv,col(i,:));
    pause(.01)
    hold on
    fileNameOBJ = sprintf('%sdSamp%d_%s_%d.obj',objDir,downSamp,tag,obName);
    fileNameSTL = sprintf('%sdSamp%d_%d.stl',objDir,downSamp,obName);
    %STLWRITE(FILE, FACES, VERTICES)
    stlwrite(fileNameSTL,fv.faces,fv.vertices);
    vertface2obj(fv.vertices,fv.faces,fileNameOBJ,obName);
    toc
    %     cellDat(i).subs = sub;
    %     cellDat(i).fv = fv;
    end
%    disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList{i},i,length(cellList)));
end
camlight
hold off










