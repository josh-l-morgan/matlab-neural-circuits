% clear all
% load('MPN.mat')
% 
% load([MPN 'obI.mat'])
% tic
% load([MPN 'dsObj.mat'])
% toc 

region = 2; %region 1 is shaft, region 2 is targeted
downSamp = 2;

seedList = [ 125];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
allEdges = obI.nameProps.edges(:,1:2);

mot = getMotifs(obI);
rgcs = mot.cel.types.rgcs;
tcrs = mot.cel.types.tcrs;
lins = mot.cel.types.lins;

Pre = postTo(allEdges,seedList);
Post = preTo(allEdges,seedList);

if region == 1
    isLocal = [ 9078 805 399 426 436 423 9003 434 819 349 492 9193];% 221 353 402]; %% shaft synapses1
else
   isLocal = [125 9013 9014 9015 9016 156 394]; 
end
% alls0 394 9014m   (12275, 15128, 7418) ,   (11760, 15418, 7736) 9013, 394 part 9015

TCR = setdiff(intersect([Post(:,1)],tcrs),seedList)';
RGC = setdiff(intersect([Pre(:,1)],rgcs),seedList)';
LINout = setdiff(intersect([Post(:,1)],lins),seedList)';
LINin = setdiff(intersect([Pre(:,1)],lins),seedList)';

synAnc = getSynAnchors(obI,125, TCR);
synBallzTCR = ballz(synAnc,3);
synAnc = getSynAnchors(obI, RGC,125);
synBallzRGC = ballz(synAnc,3);
synAnc = getSynAnchors(obI,125, LINout);
synBallzLINout = ballz(synAnc,3);
synAnc = getSynAnchors(obI, LINin,125);
synBallzLINin = ballz(synAnc,3);

if 0;%region == 2
    synAnc = getSynAnchors(obI,156,394);
    synBalzMore = ballz(synAnc,5);
else
   synBalzMore = {[]}; 
end

linked = [];
%cellList = [125 385 ]
cellList = {125 cat(1,synBallzTCR{:}) cat(1,synBallzRGC{:}) ...
    cat(1,synBallzLINout{:}) cat(1,synBallzLINin{:}) cat(1,synBalzMore{:})};
col = [1 0 0; 0 0 1; 0 1 0; 1 .5 0; 1 .5 0; 1 1 1];
cellAlpha = [1 .7 .7 .7 .7 .7];


% 
% cellList = { cat(1,synBallz{:})};
% col = [ 1 0 0];
% cellAlpha = [ .6]



%%

renderOb = 0;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end


% manualSyn = [(11763, 15425, 7696);(12120, 15265, 7468); (12167, 15199, 7447);...
%     (12248, 15206, 7445);(12306, 15145, 7395);(12405, 15156, 7390)];


%target = [120 110 110]*8; % ( X Y Z)
%rad = 300;
%crop = [target- rad; target + rad];
%crop     = [ 660   880   680;           1100   1100  1000]; %[ z x y] ; %faciculation
%crop     = [ 760   880   580;           1100   1200  1000]; %[ z x y]
%crop     = [1800  950 1050;          1880 1050 1150]; %[ z x y
%target = [908 490 555.5]*2; % ( X Y Z)

if region == 1
    target = dsAnchors([18245  20996  3042],obI,[2 1 3]);
    crop = [target - [100 200  120]; target + [100 200 80]];
    dsAnchorsReverse(target,obI,[2 1 3])
else
    crop = [1750  900 950;          1920 1100 1250]; %[ z x y
    crop = [1650  920 1080;          1960 1050 1250]; %[ z x y
    
end
clf
% cellList = 125;
l = lightangle(0,45) ;

for i = 1:length(cellList)
    subCell = names2Subs(obI,dsObj,cellList(i));
    sub = subCell{1};
    obName = cellList{i};
    if ~isempty(sub)
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
        disp(sprintf('no points on %d',i))
    else
    fv = subVolFV(smallSub,[],renderOb);
    [p] = renderFV(fv,col(i,:),cellAlpha(i));

    
         p.DiffuseStrength = 1;
         p.AmbientStrength = .7;
    
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
end
hold off


%% movie
tag = 'testMove';
frames = 360;
el = ones(frames,1) * 30;
az = 1:360/frames:360;
obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_Clear_4\'
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
savefig([obMovDir tag '.fig'])

%
% cam2 = light
% cam3 = camlight('headlight')
% set(cam2,'Position',[1 1 1])
shouldWrite = 0;

imageName = sprintf('%sspringRun_%s%05.0f_syns_01h.png',obMovDir,tag,i);

if region == 1
    % [az el] = view
    
    view(19.8, -62.8)
    lightangle(l,19.8, -22.8)
    pause(.01)
    axis off
    
else
    
    view([82.4 32.4])
    lightangle(l,82.4,110)
    pause(.01)
    axis off
    
end
return
if shouldWrite
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
end

while 0
    for i = 94%1:frames;
        
        view([az(i) el(i)])
        lightangle(l,az(i)+10, 50)
        pause(.01)
        axis off
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
    break
end




