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

allEdges = obI.nameProps.edges;
postTarg = preTo(allEdges,125);
preTarg = postTo(allEdges,125);
cellList = [preTarg(:,1)' postTarg(:,1)' 125]
%cellList = {125 'frag 125'};
% 

group1 = setdiff(preTarg(:,1)',[125 0]);
group2 = setdiff(postTarg(:,1)',[125 0]);
% 
% group1 = group1(1);
% group2 = group2(1);

group1 = [399  ];
group2 = 426;
cellList = [125 group1 group2];
colNum = length(cellList)-1;
colMap = hsv(256);
rainBow = colMap(ceil([1:colNum] * 255/(colNum)),:);
rainBow = rainBow(randperm(size(rainBow,1)),:);
%col = [rainBow; [1 1 1]]
col = [1 0 0;repmat([0 .4 1],[length(group1) 1]);repmat([.4 0 1],[length(group2) 1])]
alph = [1;repmat([.5],[length(group1) 1]);repmat([.5],[length(group2) 1])]

%col = [1 1 1];
%%

renderOb = 1;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end



downSamp = 4;
target = dsAnchors([18245  20996  3042],obI,[2 1 3]);
%target = [908 490 555.5]*2; % ( X Y Z)
rad = 300;
crop = [target - [100 200  100]; target + [100 200 100]];
%crop = [target- rad; target + rad];
% crop     = [1600  850 850;          1950 1100 1430]; %[ z x y
% crop     = [1750  900 950;          1920 1100 1250]; %[ z x y

clf
% cellList = 125;
l = lightangle(0,45) ;

for i = 1:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
    obName = cellList(i);
    if ~isempty(sub)
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
        disp(sprintf('no points on %d',cellList(i)))
    else
    fv = subVolFV(smallSub,[],0);
    [p] = renderFV(fv,col(i,:));
    view([0 0])
    
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
    end
    disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
end

hold off


%% movie
tag = 'testMove';
frames = 360;
el = ones(frames,1) * 30;
az = 1:360/frames:360;
obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_shaft01\'
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
savefig([obMovDir tag '.fig'])

% 
% cam2 = light
% cam3 = camlight('headlight')
% set(cam2,'Position',[1 1 1])
shouldWrite = 0;



while 1
for i = 1:frames;
    
view([az(i) el(i)])
lightangle(l,az(i)+10, 50)
axis off
pause%(.01)

set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])

%runSprings(springDat,allResults{1})
set(gcf, 'InvertHardCopy', 'off');
imageName = sprintf('%sspringRun_%s%05.0f.png',obMovDir,tag,i);
%print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')

if shouldWrite
    print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
end


end
if shouldWrite, break,end
end






