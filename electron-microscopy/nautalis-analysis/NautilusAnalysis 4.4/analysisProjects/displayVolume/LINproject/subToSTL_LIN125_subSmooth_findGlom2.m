clear all

obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_Glom2_SynON_scaleBar2\'
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
% else
%     'directory already exists'
%     return
% end


load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


region =2;
flipDim = [3 2 1 ];

markersOn = 1;
shouldWrite = 0;
onlyTest = 0;

if region == 1
    downSamp = 4;
    renderProps.smooth = 2;
    renderProps.resize = 1;
    renderProps.smoothPatch = 2;
else
    downSamp = 1;
    renderProps.smooth = 1;
    renderProps.resize = 1;
    renderProps.smoothPatch = 1;
end


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

if region == 1
    isLocal = [RGC];%[ 9078 805 399 426 436 423 9003 434 819 349 492 9193];% 221 353 402]; %% shaft synapses1
else
    isLocal = [RGC TCR LINout LINin];
    %isLocal = [125 9014 9018 9019 9020 9021  156 527 ]; %[125 9014 9018 9019 9020 9021 150 156 527 397];
    isLocal = RGC;
end
group1 = intersect(RGC, isLocal);
group2 = intersect(TCR, isLocal);
group3 = intersect([LINin LINout], isLocal);
group4 = [];%[227]
%
% group1 = 156
% group2 = 527
% group3 = 9019;
%
% group1 = group1(1);
% group2 = group2(1);

cellList = [125 group1 group2 group3 group4];
%cellList = [125];
% colNum = length(cellList)-1;
% colMap = hsv(256);
% rainBow = colMap(ceil([1:colNum] * 255/(colNum)),:);
% rainBow = rainBow(randperm(size(rainBow,1)),:);
%col = [rainBow; [1 1 1]]
groupCol = [1 0 0; 0 1 0; 0 0 1; 1 .5 0; .7 .7 .7];

if region == 1
    groupAlph = [1 1 1 1 1];
else
    groupAlph = [1 .5 .5 .5 .5 .5];
    groupAlph = [1 1 1 1 1 1];
    groupAlph = [1 .6 1 1 1 1];
    
end

if region == 1
    group2Col = repmat(groupCol(2,:),[length(group1) 1])
else
    group2Col = repmat(groupCol(2,:),[length(group1) 1])



end

col = [groupCol(1,:);group2Col; ...
    repmat(groupCol(3,:),[length(group2) 1]); repmat(groupCol(4,:),[length(group3) 1]);...
    repmat(groupCol(5,:),[length(group4) 1])];
alph = [groupAlph(1);repmat(groupAlph(2),[length(group1) 1]);repmat(groupAlph(3),[length(group2) 1]);...
    repmat(groupAlph(4),[length(group3) 1]); repmat(groupAlph(5),[length(group4) 1])];

%col = col + rand(size(col))*.5;
%col(2:5,:) = col(2:5,:) + [0 0 0; .7 0 0; 0 0 .7; .7 0 .7];
col(1,:) = [1 0 0]
col(col>1) = 1;

%col = [1 0 0; .5 1 0 ;0 .8 0; 0 1 .5; .2 1 .2; .8 1 .8; .5 0 1; 0 .5 1]

%col = [1 1 1];




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
    crop = [1 1 1; max(cat(1,dsObj.subs),[],1)];
else
    target = dsAnchors([2757, 14516, 8956],obI,[2 1 3]);
    crop = [target - [100 200  120]; target + [100 200 80]];
	crop = [target - [300 300  320]; target + [300 300 300]];

    %crop =  [1758    0951    1004 ; 1858    1050    1095];
    dsAnchorsReverse(target,obI,[2 1 3])
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
        if exist('crop','var')
            useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
            
            sub = sub(useSub,:);
        end
        sub = sub(:,flipDim);
        sub = sub(:,[2 1 3]);
        if isfield(renderProps,'resize')
            sub = sub* renderProps.resize;
        end
        sub = sub/downSamp;%shrinkSub(sub(:,[2 1 3]),downSamp);
        %sub = sub(:,flipDim);
        %scatter3(sub(:,1),sub(:,2),sub(:,3),10,'o','filled','w')
        
        c = 0;
        cone = patchShape('cone',6,32);
        cone.vertices = cone.vertices(:,[2 3 1])/downSamp * .6;
        cone.vertices(:,2) = cone.vertices(:,2) * -1;
        if isfield(renderProps,'resize')
            cone.vertices = cone.vertices * renderProps.resize;
        end
        cones = cone;
        for i = 1:size(sub,1)
            shiftSA = sub(i,:);
            %shiftSA = shiftSA(:,flipDim);
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
cellsShown = {};
cellsNotShown = {};
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
        smallSub = smallSub(:,flipDim);
        if onlyTest
            smallSub = smallSub(1:400,:);
        end
        trackCells = cat(1,trackCells,[obName size(sub,1)]);
        tic
        if isempty(smallSub)
            disp(sprintf('no points on %d',cellList(i)))
            cellsNotShown = cat(1,cellsNotShown,{cellList(i)});
        else
            cellsShown = cat(1,cellsShown,{cellList(i)});

            fv = subVolFV(smallSub,[],renderProps);
            [p] = renderFV(fv,col(i,:),alph(i));
            view([0 0])
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
    end
    disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
end
%aL = lightangle(0,45) ;

trackCells = trackCells(trackCells(:,2)>0,:)

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
    
    startAngle =[20 32.4];
    startAngle =[20 -22.4];
    view(startAngle)
    lightangle(aL,19.8, -22.8)
    pause(.01)
    axis off
    
else
    startAngle = [82.4 32.4];
    startAngle = [68.2 24.4];
    startAngle = [75.8 8.4];
    startAngle = [-163 58];
    startAngle = [ 61    1]

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


%%%

if 1 %scale bar
    
    
    scaleLength = 10;
    xyzC = mean(crop,1);
    xyz1 = mean(crop,1);
    xyz2 = mean(crop,1);
    xyz2(2) = xyz2(2)+ scaleLength/obI.em.dsRes(1);
    xyz1 = shrinkSub(xyz1,downSamp);
    xyz2 = shrinkSub(xyz2,downSamp);
    xyzC = shrinkSub(xyzC,downSamp);
    
    if region == 2
        
        
        xyz1(1) = xyz1(1) + 50;
        xyz1(2) = xyz1(2) + -10;
        xyz1(3) = xyz1(3) + 0;
        xyz2(1) = xyz2(1) + 50;
        xyz2(2) = xyz2(2) + -10;
        xyz2(3) = xyz2(3) + 0;
        
        
    else
        xyz1(1) = xyz1(1) + 0;
        xyz1(2) = xyz1(2) + 50;
        xyz1(3) = xyz1(3) + 70;
        xyz2(1) = xyz2(1) + 0;
        xyz2(2) = xyz2(2) + 50;
        xyz2(3) = xyz2(3) + 70;
        
    end
    if isfield(renderProps,'resize')
        xyz1 = xyz1* renderProps.resize;
        xyz2 = xyz2* renderProps.resize;
        xyzC = xyzC* renderProps.resize;
    end
    
    
    xyzC = xyzC(:,flipDim);
    xyz1 = xyz1(:,flipDim);
    xyz2 = xyz2(:,flipDim);
    
    scaleBar = plot3([xyz1(2); xyz2(2)],[xyz1(1); xyz2(1)],[xyz1(3); xyz2(3)],'linewidth',10,'color','w')
    pause(1)
    delete(scaleBar)
    dot = scatter3(xyzC(2),xyzC(1),xyzC(3),'filled','w')
end

%%%
for i = 1:length(az);
    i
    delete(scaleBar)
    
    view([az(i)+startAngle(1) startAngle(2)])
    view([az(i)+startAngle(1) startAngle(2)])
    
    lightangle(aL,az(i)+10+startAngle(1),startAngle(2)+30)
    
    %%rotate scalebar
    
    mz = makehgtform('zrotate',(az(i)+startAngle(1))*2*pi/360);
    %mx = makehgtform('yrotate',0);
    m = mz;
    
    
    xyz1r = [xyz1-xyzC 1];
    xyz1r = xyz1r * m;
    xyz2r = [xyz2-xyzC 1];
    xyz2r = xyz2r * m;
    xyz1r = xyz1r(1:3)+xyzC;
    xyz2r = xyz2r(1:3) + xyzC;
        
    scaleBar = plot3([xyz1r(2); xyz2r(2)],[xyz1r(1); xyz2r(1)],[xyz1r(3); xyz2r(3)],'linewidth',10,'color','w');
        
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




