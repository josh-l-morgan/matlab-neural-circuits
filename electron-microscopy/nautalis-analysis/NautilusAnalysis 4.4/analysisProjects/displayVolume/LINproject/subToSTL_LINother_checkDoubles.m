clear all
load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

cellList = [unique([9017  527 9018 156 527 397 9013 394  527 394 ...
    9014 9015 156 9015 9016 395 396 9019 9020 9021]) 125];

cellList = [ 156 527 125];

[postTarg] = getList_125rand;
postAx = [447:456 10034:10036 160 303 184];
cellList = [125 postAx postTarg'];

cellList = {125 9026};

col = [ 1 0 0 ; 0 1 1];
alph = [.1 1];


renderOb = 0;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end



downSamp = 4;

target = [908 490 555.5]*2; % ( Y X Z)
rad = 300;
%crop = [target- rad; target + rad];
%crop     = [2000 400 600;          2500 600 700];
clf
for i = 1:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
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
        disp(sprintf('no points on %d',cellList{i}))
    else
    fv = subVolFV(smallSub,[],renderOb);
    renderFV(fv,col(i,:),alph(i));
    pause(.01)
    hold on
    fileNameOBJ = sprintf('%sdSamp%d_%s_%d.obj',objDir,downSamp,tag,obName);
    fileNameSTL = sprintf('%sdSamp%d_%d.stl',objDir,downSamp,obName);
    %STLWRITE(FILE, FACES, VERTICES)
    %stlwrite(fileNameSTL,fv.faces,fv.vertices);
    vertface2obj(fv.vertices,fv.faces,fileNameOBJ,obName);
    toc
    %     cellDat(i).subs = sub;
    %     cellDat(i).fv = fv;
    end
    disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList{i},i,length(cellList)));
end
camlight
hold off










