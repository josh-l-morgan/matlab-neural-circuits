



if 0
    clear all
    
    obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_TargFac\'
    if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
    % else
    %     'directory already exists'
    %     return
    % end
    
    load('MPN.mat')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
    
    sm = addDatToSynMat(obI)
    sm  = getTopoEucDistBetweenSyn(sm);
    sm = getTopoEucDistBetweenSkelAndSyn(sm);
    sm = labelShaftSkel(sm);
end

region = 2;
downSamp = 1;
markersOn = 0;
shouldWrite = 0;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 0;
onlyTest = 0;



%% Convert pos [Y X Z]
targL{1} = [834.1 166.4  606.8];
targL{2} = [1057 868.9 198];


%Top
targL{3} = [391.9 395 417.4];
targL{4} = [434.7 389 345]; % same as 3
targL{5} = [385.6 474.1 429.1];
targL{6} = [499.1 511.5 391.1];
targL{7} = [554.9 467.3 414];
targL{8} = [374.2 354.2 228.9];
targL{9} = [462.2 411.8 244.7]; % same as 8
targL{10} = [458.9 461 266.2]; % same as 8
targL{11} = [550.7 380.2 281.1]; %same as 8
targL{12} = [546 448.4 289.1]; % same as 8
targL{13} = [557.7 465.9 406.7]; 
targL{14} = [731.3 450 338.3];
targL{15} = [799 467.8 289.6];

targL{16} = [1091 918.9 161.7];
targL{17} = [1138 844.1 293.7];
targL{18} = [1109 924.1 313.2];
targL{19} = [1156 958.9 457.5];

targL{20} = [1198 581 117.4]; % might be too small
targL{21} = [ 1233 612.4 209.1]; % same as 20
targL{22} = [1140 525.3 35.07]; % same as 20

targL{23} = [911.9 506.5 526];

targL{24} = [927.3 472 596.6];
targL{25} = [975.1 496.1 647.7];

targL{26} = [1036 447.1 559.7];
targL{27} = [1053 310.1 672.4];
targL{28} = [1044 359 613.7];

targL{29} = [835.2 166.8 605.7]; %Main rendered targeted process
targL{30} = [855.4 74.19 671.1];
targL{31} = [931.8 125.7 622.2];
targL{32} = [945.4 187 606];
targL{33} = [909.6 96.06 656.9];

targL{35} = [1178 292 430.2];
targL{36} = [1182 252.6 557];
targL{37} = [1145 248.1 611.4];
targL{38} = [1058 240.6 626.1];
targL{39} = [1174 182.7 649.1];
targL{40} = [1195 120 688.6];

targL{41} = [1266 165.1 595.2];
targL{42} = [1253 219.9 609.4];
targL{43} = [1302 229.1 628.2];
targL{44} = [1278 189.2 651.6];
targL{45} = [1328 144.1 686.1];

targL{46} = [ 1159 571.2 199.6];

if 1
   
    for i = 1:length(targL)
        i
        hold on
        set(gca,'Clipping','off')

        markTarg = targL{i};
        if length(markTarg) == 3
            scatter3(markTarg(2),markTarg(1),markTarg(3),'o','w','filled')
        else
            markTarg
            i
        end
        %pause
    end

    view([00 0])
end


for t = 1:length(targL);
    
target = targL{t};
YXZdownSamp = 2;

targUM = target * .2 * YXZdownSamp;
targVox = target * YXZdownSamp;
viewSize = [100 100 100];
crop = [targVox - viewSize; targVox + viewSize];
cropUM = crop * 0.2;

%allPos2 = allPos;
% allPos2 = [Y X Z] * YXZdownSamp;
%     target = dsAnchors([18245  20996  3042],obI,[2 1 3]);
%     dsAnchorsReverse(target,obI,[2 1 3])


%% Synapse filter

targCell = 125;
filtPre = sm.pre == targCell;
filtPost = sm.postClass == 2;
filtPos = ((sm.pos(:,1)>cropUM(1,1)) & (sm.pos(:,1)<cropUM(2,1)) & ...
    (sm.pos(:,2)>cropUM(1,2)) & (sm.pos(:,2)<cropUM(2,2)) & ...
    (sm.pos(:,3)>cropUM(1,3)) & (sm.pos(:,3)<cropUM(2,3)));
isSyn = find(filtPre & filtPost & filtPos);
isSyn = find(filtPre  & filtPos);



cellList = unique([sm.post(isSyn) sm.pre(isSyn)]);



%%

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
    isLocal = [ 9078 805 399 426 436 423 9003 434 819 349 492 9193];% 221 353 402]; %% shaft synapses1
else
    isLocal = [125 9013 9014 9015 9016 156 394];
end
group1 = intersect(RGC, cellList);
group2 = intersect(TCR, cellList);
group3 = intersect([LINin LINout], cellList);
group4 = [];%[227]
%
% group1 = 156
% group2 = 527
% group3 = 9019;
%
% group1 = group1(1);
% group2 = group2(1);

cellList = [125 group1(:)' group2(:)' group3(:)' group4(:)'];
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
    group2Col = [0 1 0; 0 1 .8; .8 1 0; .4 1 0];
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
        if onlyTest
            smallSub = smallSub(1:400,:);
        end
        trackCells = cat(1,trackCells,[obName size(sub,1)]);
        tic
        if isempty(smallSub)
            disp(sprintf('no points on %d',cellList(i)))
        else
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
az = 0; %0:360/frames:359;
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


for i = 1:length(az);
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


end


