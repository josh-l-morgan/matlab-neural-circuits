clear all
load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

seedList = [ 125];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
allEdges = obI.nameProps.edges(:,1:2);

Pre = postTo(allEdges,seedList);
Post = preTo(allEdges,seedList);
TCR = setdiff([Post(:,1)],seedList)';
RGC = setdiff([Pre(:,1)],seedList)';


synAnc = getSynAnchors(obI,125, 385);
synBallz = ballz(synAnc,12);




linked = [];
%cellList = [125 385 ]
cellList = {125 385 'frag 385 branch 125'  'frag 125 branch 385' cat(1,synBallz{:})};
col = [0 .8 0; .8 0 .8; 1 0 1; 0 1 0; 1 0 0];
cellAlpha = [.3 .3 1 1 .4]

cellList = {125 385 'frag 385 branch 125'  'frag 125 branch 385' cat(1,synBallz{:})};
col = [0 .8 0; .8 0 .8; 1 0 1; 0 1 0; 1 0 0];
cellAlpha = [.3 .3 1 1 .5]

% 
% cellList = { cat(1,synBallz{:})};
% col = [ 1 0 0];
% cellAlpha = [ .6]



%%

renderOb = 0;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end



downSamp = 4;

target = [120 110 110]*8; % ( X Y Z)
rad = 300;
%crop = [target- rad; target + rad];
%crop     = [ 660   880   680;           1100   1100  1000]; %[ z x y] ; %faciculation
crop     = [ 760   880   580;           1100   1200  1000]; %[ z x y]

clf
% cellList = 125;
l = lightangle(0,45) ;

for i = 1:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
    obName = cellList{i};
    
    if ~ischar(obName)
        if length(obName>1)
            obName = sprintf('object_%d',i);
        end
    end
    
    if iscell(obName); obName = obName{1};end
    if exist('crop','var')
        useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
            (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
            (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
        sub = sub(useSub,:);
        
    end
    smallSub = shrinkSub(sub,downSamp);
    smallSub = smallSub(:,[2 3 1]);

    tic
    if isempty(smallSub)
        disp(sprintf('no points on %d',cellList(i)))
    else
    fv = subVolFV(smallSub,[],renderOb);
    [p] = renderFV(fv,col(i,:),cellAlpha(i));

    view([0 0])
set(gca,'Zdir','reverse')
    axis off
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
    %disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
end

hold off


%% movie
tag = 'testMov2e2';
frames = 180;
el = ones(frames,1) * 0;
az = 1:360/frames:360;
obMovDir = 'D:\LGNs1\Analysis\movies\seed125_11\'
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
savefig([obMovDir tag '.fig'])

% 
% cam2 = light
% cam3 = camlight('headlight')
% set(cam2,'Position',[1 1 1])
shouldWrite = 1;



while 0
for i = 1:frames;
    
view([az(i) el(i)])
axis off
set(gca,'Zdir','reverse')
lightangle(l,az(i)-20, -30)
pause(.01)

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






