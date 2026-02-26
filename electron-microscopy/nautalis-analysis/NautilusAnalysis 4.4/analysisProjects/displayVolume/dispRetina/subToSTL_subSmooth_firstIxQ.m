clear all

load('MPN.mat')
load('WPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


obMovDir = [WPN '\movies\IxQ\firstRenderVG3\'];
if ~exist(obMovDir,'dir'),mkdir(obMovDir),end


region =1;
flipDim = [1 3 2 ];

markersOn = 1;
synSegOn = 1;
shouldWrite = 0;
onlyTest = 0;
rotate = 0;
showScale = 0;
fullContext = 0;
randCol = 0;

if region == 1
    downSamp = 1;
    renderProps.smooth = 0;
    renderProps.resize = 2;
    renderProps.smoothPatch = 10;
else
    downSamp = 1;
    renderProps.smooth = 4;
    renderProps.resize = 2;
    renderProps.smoothPatch = 2;
    %renderProps.dilate = 1;
end

%% Get cell info
allEdges = obI.nameProps.edges;
allIds = unique(obI.nameProps.cellNum);
allSegs = 1:length(dsObj);
%allCells = obI.cell.name;


%% Pick cells
% mot = getMotifs(obI);
% allCells = mot.cel.cells;
% rgcs = mot.cel.types.rgcs;
% tcrs = mot.cel.types.tcrs;
% lins = mot.cel.types.lins;
%
seedList = 4;
Post = preTo(allEdges,seedList);
Pre = postTo(allEdges,seedList);
%
% TCR = setdiff(intersect([Post(:,1)],tcrs),seedList)';
% RGC = setdiff(intersect([Pre(:,1)],rgcs),seedList)';
% LINout = setdiff(intersect([Post(:,1)],lins),seedList)';
% LINin = setdiff(intersect([Pre(:,1)],lins),seedList)';

if region == 1
    isLocal = [allIds ];%[9079 353 468 221 362 5119];%[ 9078 805 399 426 436 423 9003 434 819 349 492 9193];% 221 353 402]; %% shaft synapses1
    group1 = [];
    group2 = [];
    group3 = seedList;
    group4 = setdiff(Pre(:,1),[seedList 0])';%
    group5 = setdiff(Post(:,1),[seedList 0])';
    group6 = setdiff(allIds,[group1 group2 group3 group4 group5]);
else
    isLocal = [RGC TCR LINout LINin];
    %isLocal = [125 9014 9018 9019 9020 9021  156 527 ]; %[125 9014 9018 9019 9020 9021 150 156 527 397];
    isLocal = RGC;
    group1 = {'rgcD' 'rgcF' 'rgcG' 'rgcJ' 'rgcL' 'rgcP' 'rgcQ' 'rgcX' 'rgcAE'};%intersect(RGC, isLocal);
    group2 = {'tcrB' 'tcrC' 'tcrI' 'tcrO' 'tcrY'};%intersect(TCR, isLocal);
    group3 = {'linH' 'linM' 'linR' 'linT' 'linZ' 'linAF' 'linAG'};%intersect([LINin LINout], isLocal);
    group4 = {'unkK' 'unkS' 'unkU' 'unkV' 'unkAA' 'unkAB' 'unkAC' 'unkAD' };%[227]
    group5 = cat(2,num2cell(setdiff(allCells,125)),{'tcr C'})
    group6 = [];
    
end


cellList = cat(2,group1, group2, group3, group4, group5, group6);

if randCol
    colNum = length(cellList);
    colMap = hsv(256);
    rainBow = colMap(ceil([1:colNum] * 255/(colNum)),:);
    rainBow = rainBow(randperm(size(rainBow,1)),:);
    groupCol = [rainBow];
    groupAlph = [1 1 1 1 1 1];
    
else
    groupCol = [1 1 1 ; 1 1 1; 1 1 1; 0 1 0; 0 0 1; .7 .7 .7; .4 .4 .4; 0 .5 1];
    groupAlph = [.2 .2 .4 .2 .2 .1];
end



col = [repmat(groupCol(1,:),[length(group1) 1]);repmat(groupCol(2,:),[length(group2) 1]);
    repmat(groupCol(3,:),[length(group3) 1]); repmat(groupCol(4,:),[length(group4) 1]);...
    repmat(groupCol(5,:),[length(group5) 1]);repmat(groupCol(6,:),[length(group6) 1])];
alph = [repmat(groupAlph(1),[length(group1) 1]);repmat(groupAlph(2),[length(group2) 1]);repmat(groupAlph(3),[length(group3) 1]);...
    repmat(groupAlph(4),[length(group4) 1]); repmat(groupAlph(5),[length(group5) 1]); repmat(groupAlph(6),[length(group6) 1])];

col(col>1) = 1;




%%
renderOb = 1;
tag = 'testCrop';
objDir = [WPN 'stlFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end

if region == 1
    target = dsAnchors([14938, 20514, 2869],obI,[2 1 3])
    crop = [1 1 1; max(cat(1,dsObj.subs),[],1)]*8;
    %crop = [target - [40 50  50]; target + [50 80 40]];
    %crop = [2361 995 310; 2562 1395 5103];
else
    target = dsAnchors([2757, 14516, 8956],obI,[2 1 3]);
    crop = [target - [100 200  120]; target + [100 200 80]];
    dsAnchorsReverse(target,obI,[2 1 3])
end
%dsAnchorsReverse(crop,obI,[2 1 3])
%datPos = dsAnchorsReverse(mean(crop,1),obI,[2 1 3])

clf

%% Draw Markers

[marker(1).sub marker(1).segIDs] = getSynAnchors(obI, allIds, group3, [2 3]);
[marker(2).sub marker(2).segIDs] = getSynAnchors(obI, allIds, group3, 1);
[marker(3).sub marker(3).segIDs] = getSynAnchors(obI, group3, allIds );
[marker(4).sub marker(4).segIDs] = getSynAnchors(obI, group3, allIds, [2 3]);

if markersOn
    
    markerCol = [0 1 0; 1 0 0; .2 .2 1; 0 1 1];
    markerAlpha = [.3 .3 .3 .3 .3 .3];
    
    for m = 1:length(marker);
        sub = marker(m).sub;
        if ~isempty(sub)
            if exist('crop','var')
                useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                    (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                    (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
                sub = sub(useSub,:);
            end
            sub = sub(:,flipDim);
            sub = sub(:,[2 3 1]);
            if isfield(renderProps,'resize')
                sub = sub* renderProps.resize;
            end
            sub = sub/downSamp;%shrinkSub(sub(:,[2 1 3]),downSamp);
            
            c = 0;
            cone = patchShape('cone',6*2,16*2);
            cone.vertices = cone.vertices(:,[2 3 1])/downSamp * .6;
            cone.vertices(:,1) = cone.vertices(:,1) * -1;
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
                cs.faceAlpha = markerAlpha(m);
                patch(cs)
                hold on
            end
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
            if (fullContext == 0) | (i>1)
                useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                    (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                    (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
                sub = sub(useSub,:);
            end
        end
        smallSub = shrinkSub(sub,downSamp);
        smallSub = smallSub(:,flipDim);
        if onlyTest
            smallSub = smallSub(1:400,:);
        end
        %trackCells = cat(1,trackCells,[obName size(sub,1)]);
        tic
        if isempty(smallSub)
            %disp(sprintf('no points on %d',cellList(i)))
            cellsNotShown = cat(1,cellsNotShown,{cellList(i)});
        else
            cellsShown = cat(1,cellsShown,{cellList(i)});
            
            fv = subVolFV(smallSub,[],renderProps);
            [p] = renderFV(fv,col(i,:),alph(i));
            view([0 0])
            axis off
            pause(.01)
            hold on
%             fileNameOBJ = sprintf('%sdSamp%d_%s_%d.obj',objDir,downSamp,tag,obName);
%             fileNameSTL = sprintf('%sdSamp%d_%d.stl',objDir,downSamp,obName);
            %STLWRITE(FILE, FACES, VERTICES)
            %stlwrite(fileNameSTL,fv.faces,fv.vertices);
            %vertface2obj(fv.vertices,fv.faces,fileNameOBJ,obName);
            toc
            %     cellDat(i).subs = sub;
            %     cellDat(i).fv = fv;
        end
    end
    % disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
end
%aL = lightangle(0,45) ;

%trackCells = trackCells(trackCells(:,2)>0,:)



%% Draw Synapses

synSegCol =  markerCol ;
synSegAlph  = [1 1 1 1 1 1 1 1  ];

if synSegOn
    for i = 1:length(marker)
        
        sub = cat(1,dsObj(marker(i).segIDs).subs);
        
        if ~isempty(sub)
            if iscell(obName); obName = obName{1};end
            if exist('crop','var')
                if (fullContext == 0) | (i>1)
                    useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                        (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                        (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
                    sub = sub(useSub,:);
                end
            end
            smallSub = shrinkSub(sub,downSamp);
            smallSub = smallSub(:,flipDim);
            if onlyTest
                smallSub = smallSub(1:400,:);
            end
            %trackCells = cat(1,trackCells,[obName size(sub,1)]);
            tic
            if isempty(smallSub)
                %disp(sprintf('no points on %d',cellList(i)))
                cellsNotShown = cat(1,cellsNotShown,{cellList(i)});
            else
                
                
                fv = subVolFV(smallSub,[],renderProps);
                [p] = renderFV(fv,synSegCol(i,:),synSegAlph(i));
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
        % disp(sprintf('finished rendering cell %d.  (%d of %d)',cellList(i),i,length(cellList)));
    end
end

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
    startAngle =[0 0];
    view(startAngle)
    lightangle(aL,19.8, -22.8)
    pause(.01)
    axis off
    
else
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

if showScale %scale bar
    
    scaleLength = 10;
    xyzC = mean(crop,1);
    xyz1 = mean(crop,1);
    xyz2 = mean(crop,1);
    xyz2(2) = xyz2(2)+ scaleLength/obI.em.dsRes(1);
    
    if 0 % define points
        a = [14794, 20455, 2869]
        b = [15034, 20619, 2869]
        %
        %     a = [14794, 20619, 2869]
        %     b = [15034, 20435, 2869]
        
        xyz1 = dsAnchors(a,obI,[2 1 3]);
        xyz2 = dsAnchors(b,obI,[2 1 3]);
    end
    
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
       
        
    end
    if isfield(renderProps,'resize')
        xyz1 = xyz1* renderProps.resize;
        xyz2 = xyz2* renderProps.resize;
        xyzC = xyzC* renderProps.resize;
    end
    
    
    xyzC = xyzC(:,flipDim);
    xyz1 = xyz1(:,flipDim);
    xyz2 = xyz2(:,flipDim);
    
    
    scaleBar = plot3([xyz1(2); xyz2(2)],[xyz1(1); xyz2(1)],[xyz1(3); xyz2(3)],'linewidth',2,'color','w')
    dot = scatter3(xyzC(2),xyzC(1),xyzC(3),'filled','w')
    
    delete(scaleBar)
    delete(dot)
end

%%%

if ~rotate
    view([0 0])
    lightangle(aL,90,90)
else
    
    %%
    for i = 1:length(az);
        i
        if exist('scaleBar'), delete(scaleBar); end
        
        view([az(i)+startAngle(1) startAngle(2)])
        lightangle(aL,az(i)+10+startAngle(1),startAngle(2)+30)
        
        %%rotate scalebar
        mz = makehgtform('zrotate',(az(i)+startAngle(1))*2*pi/360);
        %mx = makehgtform('yrotate',0);
        m = mz;
        
        if showScale
            xyz1r = [xyz1-xyzC 1];
            xyz1r = xyz1r * m;
            xyz2r = [xyz2-xyzC 1];
            xyz2r = xyz2r * m;
            xyz1r = xyz1r(1:3)+xyzC;
            xyz2r = xyz2r(1:3) + xyzC;
            
            
            scaleBar = plot3([xyz1r(2); xyz2r(2)],[xyz1r(1); xyz2r(1)],[xyz1r(3); xyz2r(3)],'linewidth',10,'color','w');
        end
        
        pause(.01)
        axis off
        set(gcf,'PaperUnits','points','PaperPosition',[1 1 512 512])
        %runSprings(springDat,allResults{1})
        set(gcf, 'InvertHardCopy', 'off');
        imageName = sprintf('%srot_%05.0f.png',obMovDir,i);
        if shouldWrite
            print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
        end
        
    end
end


