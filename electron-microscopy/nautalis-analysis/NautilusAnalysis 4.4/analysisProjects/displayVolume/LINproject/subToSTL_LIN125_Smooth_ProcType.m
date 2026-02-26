
obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_NeuriteType\'
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
% else
%     'directory already exists'
%     return
% end

downSamp = 2;


%% Get shaft

MPNshaft = 'D:\LGNs1\Export\export_joshm_dendType_cE\'
load([MPNshaft 'obI.mat'])
load([MPNshaft 'dsObj.mat'])

checkExp = 'shaft';
isShaft = [];
shaftSub = [];
allSub = [];
axSub = [];
targSub = [];

for i = 1:length(obI.nameProps.names)
    sub = dsObj(i).subs;
    allSub = cat(1,allSub,sub);
    if sum(regexp(obI.nameProps.names{i}, 'shaft'))
        isShaft = [isShaft i];
        shaftSub = cat(1,shaftSub,sub);
    elseif sum(regexp(obI.nameProps.names{i}, 'axon'))
        axSub = cat(1,axSub,sub);
    elseif sum(regexp(obI.nameProps.names{i}, 'cellBody'))
        bodySub = cat(1,axSub,sub);
    else 
        targSub = cat(1,targSub,sub);
    end
end


herSamp = downSamp * 2;
shaftSub = double(shrinkSub(shaftSub,herSamp));
axSub = double(shrinkSub(axSub,herSamp));
targSub = double(shrinkSub(targSub,herSamp));
bodySub = double(shrinkSub(bodySub,herSamp));
allSub = double(shrinkSub(allSub,herSamp));

% 
% s(1).sub = shaftSub;
% s(2).sub = axSub;
% s(3).sub = targSub;
% s(4).sub = bodySub;

a = max(cat(1,shaftSub,axSub,targSub),[],1)

%% Get cell


load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

targCell = 125;
subCell = names2Subs(obI,dsObj,targCell);
sub = double(subCell{1});
sub = double(shrinkSub(sub,downSamp));

b = max(sub,[],1);
a./b
%s(5).sub = sub;

distThresh = 3;


s(1).sub = [];
s(2).sub = [];
s(3).sub = [];
s(4).sub = [];
s(5).sub = [];


for i = 1:size(sub,1)
    
    distShaft = min(sqrt((shaftSub(:,1)-sub(i,1)).^2 + (shaftSub(:,2)-sub(i,2)).^2 + ...
        (shaftSub(:,3)-sub(i,3)).^2));
    distAx = min(sqrt((axSub(:,1)-sub(i,1)).^2 + (axSub(:,2)-sub(i,2)).^2 + ...
        (axSub(:,3)-sub(i,3)).^2));
    distBody = min(sqrt((bodySub(:,1)-sub(i,1)).^2 + (bodySub(:,2)-sub(i,2)).^2 + ...
        (bodySub(:,3)-sub(i,3)).^2));
    
    
    if distBody <= distThresh
        s(4).sub = cat(1,s(4).sub, sub(i,:));
    elseif distAx <= distThresh
        s(2).sub = cat(1,s(2).sub, sub(i,:));
    elseif distShaft <= distThresh
        s(3).sub = cat(1,s(3).sub, sub(i,:));
    else
        s(1).sub = cat(1,s(1).sub, sub(i,:));
    end
    
       
end

s(1).sub = sub;


%%



region = 1;
markersOn = 0;
shouldWrite = 1;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 1;
onlyTest = 0;

%
% cellList = [unique([9078  527 9018 156 527 397 9013 394  527 394 ...
%     9014 9015 156 9015 9016 395 396 9019 9020 9021]) 125];
%
% cellList = [ 156 527 125];

[postTarg] = getList_125rand;
postAx = [447:456 10034:10036 160 303 184];

allEdges = obI.nameProps.edges;

mot = getMotifs(obI);
rgcs = mot.cel.types.rgcs;
tcrs = mot.cel.types.tcrs;
lins = mot.cel.types.lins;

seedList = 125;
Pre = postTo(allEdges,seedList);
Post = preTo(allEdges,seedList);

TCR = setdiff(intersect([Post(:,1)],tcrs),seedList)';
RGC = setdiff(intersect([Pre(:,1)],rgcs),seedList)';
LINout = setdiff(intersect([Post(:,1)],lins),seedList)';
LINin = setdiff(intersect([Pre(:,1)],lins),seedList)';



col = [0 1 0; 1 0 0; 0 0 1; 1 1 0; 0 1  1; 1 1 0];
alph = [1 1 1 1 1 1 1];


%%
renderOb = 1;
tag = 'testCrop';
objDir = [MPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end

% target = [908 490 555.5]*2; % ( X Y Z)
% rad = 200;
% crop = [target- rad; target + rad];
% crop     = [1600  850 850;          1950 1100 1430]; %[ z x y
% crop     = [1750  900 950;          1920 1100 1250]; %[ z x y
% crop     = [1750  900 950;          1920 1100 1250]; %[ z x y


if region == 1
%     target = dsAnchors([18245  20996  3042],obI,[2 1 3]);
%     crop = [target - [100 200  120]; target + [100 200 80]];
%     dsAnchorsReverse(target,obI,[2 1 3])
else
    crop = [1750  900 950;          1920 1100 1250]; %[ z x y
    crop = [1650  900 1080;          1960 1025 1250]; %[ z x y
    
end

%dsAnchorsReverse(crop,obI,[2 1 3])
%datPos = dsAnchorsReverse(mean(crop,1),obI,[2 1 3])


clf


%% Draw Markers
if markersOn
    groupL = [group3 125];
    % marker(1).sub = getSynAnchors(obI, group1,125);
    % marker(2).sub = getSynAnchors(obI,125, group2);
    % marker(3).sub = getSynAnchors(obI,group1, group2);
    %
    marker(1).sub = getSynAnchors(obI, group1,groupL);
    marker(2).sub = getSynAnchors(obI,groupL, group2);
    marker(3).sub = getSynAnchors(obI,group1, group2);
    marker(4).sub = getSynAnchors(obI,groupL, groupL);
    
    
    markerCol = [1 1 0; 1 0 1; 0 1 1; 1 0 0];
    %marker(4).sub = cat(1,getSynAnchors(obI,125, group3),getSynAnchors(obI,125, group3));
    
    
    for m = 1:length(marker);
        sub = marker(m).sub;
        %synAnc = dsAnchors(synAnc,obI,[2 1 3]);
        useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
            (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
            (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
       sub = sub(useSub,:);
        sub = sub(:,[2 1 3]);
        sub = sub/downSamp;%shrinkSub(sub(:,[2 1 3]),downSamp);
        
        
        scatter3(sub(:,1),sub(:,2),sub(:,3),10,'o','filled','w')
        
        c = 0;
        cone = patchShape('cone',6,40);
        cone.vertices = cone.vertices(:,[3 2 1])/downSamp;
        cone.vertices(:,2) = cone.vertices(:,2) * -1;
        cones = cone;
        for i = 1:size(sub,1)
            shiftSA = sub(i,:);
            cones(i) = cone;
            cones(i).vertices = cone.vertices + repmat(shiftSA,[size(cone.vertices,1) 1]);
            cs = cones(i);
            cs.faceColor = markerCol(m,:);
            cs.faceAlpha = 0.75;
            patch(cs)
            hold on
        end
    end
end


%% Draw cells
% cellList = 125;
aL = lightangle(0,45) ;
trackCells = [];
for i = 1:length(s)
%     subCell = names2Subs(obI,dsObj,cellList(i));
%     sub = subCell{1};
    sub = s(i).sub;
    %obName = cellList(i);
    if ~isempty(sub)
        if exist('crop','var')
            useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
            sub = sub(useSub,:);
            
        end
        smallSub = sub;%shrinkSub(sub,downSamp);
        if onlyTest
            smallSub = smallSub(1:400,:);
        end
        tic
        if isempty(smallSub)
            disp(sprintf('no points on %d',cellList(i)))
        else
            fv = subVolFV(smallSub,[],renderProps);
            [p] = renderFV(fv,col(i,:),alph(i));
            typeFV.type(i).fv = fv;
            typeFV.type(i).p = p;
            typeFV.col = col; 
            view([0 0])
            axis off
            pause(.01)
            hold on
            obName = sprintf('%04.0f',i);
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
    %disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
end
aL = lightangle(0,45) ;

save('D:\LGNs1\mergeSeg_mat\FVFiles\LIN125_procType.mat','typeFV')
%trackCells = trackCells(trackCells(:,2)>0,:)

%% movie
tag = 'testMove';
frames = 360;
el = zeros(frames,1);
az = 0:360/frames:359;
savefig([obMovDir tag '.fig'])

%
% cam2 = light
% cam3 = camlight('headlight')
% set(cam2,'Position',[1 1 1])


imageName = sprintf('%sspringRun_%s%05.0f_01h.png',obMovDir,tag,i);

if region == 1
    % [az el] = view
    startAngle = [20 22];
    startAngle = [20 32.4];
    startAngle = [15.2 -63.6];
    
    view(startAngle)
    lightangle(aL,19.8, -22.8)
    pause(.01)
    axis off
    
else
    startAngle = [82.4 32.4];
    startAngle = [68.2 24.4];
    startAngle = [75.8 8.4];
    view(startAngle)
    lightangle(aL,82.4,110)
    pause(.01)
    axis off
    %az = 0;

end


if 0
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
end


%%

startAngle = [ 0 -90 ];

for i = 1;%:length(az);
    i
    view([az(i)+startAngle(1) startAngle(2)])
    lightangle(aL,az(i)+10+startAngle(1),startAngle(2)+30)
    %lightangle(aL,82.4,110)
    pause(.01)
    axis off
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 512 512])
    %runSprings(springDat,allResults{1})
    set(gcf, 'InvertHardCopy', 'off');
    imageName = sprintf('%srot_%05.0f.png',obMovDir,i);
    %print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
    %print(gcf, imageName, '-dpng','-opengl','-r72')
    if shouldWrite
       print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
    end
    
end





