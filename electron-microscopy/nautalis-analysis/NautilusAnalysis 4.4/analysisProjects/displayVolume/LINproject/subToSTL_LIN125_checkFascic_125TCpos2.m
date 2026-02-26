
    

if 0
    clear all
    
    obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_TargFac01\'
    if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
    % else
    %     'directory already exists'
    %     return
    % end
    
    load('MPN.mat')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
    
    sm = addDatToSynMat(obI)
    sm = getTopoEucDistBetweenSyn(sm);
    sm = getTopoEucDistBetweenSkelAndSyn(sm);
    sm = labelShaftSkel(sm);
    sm = labelSubTypes(sm);

end

region = 1;
downSamp = 1;
markersOn = 1;
shouldWrite = 8;
renderProps.dilate = 0;
renderProps.smooth = 0; %(thickens
renderProps.resize = 1;
renderProps.smoothPatch = 0;
onlyTest = 0;

obName = 'heck'

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


%% Get list to extend

%%
targCell = 125;
% 
% YXZdownSamp = 2;
% % 
% % targUM = target * .2 * YXZdownSamp;
% % targVox = target * YXZdownSamp;
% viewSize = [100 100 100];
% crop = [targVox - viewSize; targVox + viewSize];
% cropUM = crop * 0.2;
% 
% 

filtPre = sm.pre == targCell;
filtPost = sm.postClass == 2;
% filtPos = ((sm.pos(:,1)>cropUM(1,1)) & (sm.pos(:,1)<cropUM(2,1)) & ...
%     (sm.pos(:,2)>cropUM(1,2)) & (sm.pos(:,2)<cropUM(2,2)) & ...
%     (sm.pos(:,3)>cropUM(1,3)) & (sm.pos(:,3)<cropUM(2,3)));
% isSyn = find(filtPre & filtPost & filtPos);
% isSyn = find(filtPre  & filtPos);
% synPos = sm.pos(isSyn,:);

% cellList = unique([sm.post(isSyn) sm.pre(isSyn)]);
checkSyn = find(filtPre & filtPost);
checkPos = sm.pos(checkSyn,:);
anchors = checkPos;
res = obI.em.res;
dSamp =  (res .* [4 4 1])./1000;
anchors(:,1) = anchors(:,1)/dSamp(1);
anchors(:,2) = anchors(:,2)/dSamp(2);
anchors(:,3) = anchors(:,3)/dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;
clear pos
pos(:,[2 1 3]) = anchors;


for t = 404%1:size(checkSyn,1);
    
    pos(t,:)
    target = sm.pos(checkSyn(t),:);
    if length(target) == 3
    YXZdownSamp = 2;
    
    targUM = target * obI.em.dsRes(1) * YXZdownSamp;
    targVox = target / obI.em.dsRes(1);
    viewSize = [60 60 60];
    crop = [targVox - viewSize; targVox + viewSize];
    cropUM = crop * obI.em.dsRes(1);
    
    %allPos2 = allPos;
    % allPos2 = [Y X Z] * YXZdownSamp;
    %     target = dsAnchors([18245  20996  3042],obI,[2 1 3]);
    %     dsAnchorsReverse(target,obI,[2 1 3])
    
    
    %% Synapse filter
    
    synPos = sm.pos(checkSyn(t),:);
    cellList = unique([sm.post(checkSyn(t)) sm.pre(checkSyn(t))]);
    
    
    %%
    
    %
    % cellList = [unique([9078  527 9018 156 527 397 9013 394  527 394 ...
    %     9014 9015 156 9015 9016 395 396 9019 9020 9021]) 125];
    %
    % cellList = [ 156 527 125];
    
    subG{1} = sm.subs(sm.subType==1,:) ;
    subG{2} = sm.subs(sm.subType==2,:) ;
    subG{3} = sm.subs(sm.subType==3,:) ;
    subG{4} = sm.subs(sm.subType==4,:);
    subCell = names2Subs(obI,dsObj,cellList(2));
    subG{5} = subCell{1};
    
    col = [0 1 0; 1 1 0; 1 0 0; 1 1 1; 0 0 1]; 
    alph = [ 1 1 1 1 1 1 1 1];
    
    %%
    renderOb = 1;
    tag = 'testCrop';
    objDir = [MPN 'stlFiles\'];
    if ~exist(objDir,'dir'),mkdir(objDir),end
    
    clf
    
    
    %% Draw Markers
    if markersOn
       
        marker(1).sub = getSynAnchors(obI, cellList(1), cellList(2));
        markerCol = [1 1 0; 1 0 1; 0 1 1; 1 0 0];
        
        
%         for m = 1:length(marker);
%             sub = marker(m).sub;
            sub = synPos ./ obI.em.dsRes(1);
            %synAnc = dsAnchors(synAnc,obI,[2 1 3]);
            useSub = ((crop(1,1)<sub(:,1)) & (crop(2,1)>sub(:,1)) & ...
                (crop(1,2)<sub(:,2)) & (crop(2,2)>sub(:,2)) & ...
                (crop(1,3)<sub(:,3)) & (crop(2,3)>sub(:,3)));
            sub = sub(useSub,:);
            sub = sub(:,[2 1 3]);
            sub = sub/downSamp * renderProps.resize ;%shrinkSub(sub(:,[2 1 3]),downSamp);
            
            scatter3(sub(:,1),sub(:,2),sub(:,3),10,'o','filled','r')
            
            c = 0;
            cone = patchShape('cone',6,40);
            cone.vertices = cone.vertices(:,[3 2 1])/downSamp*renderProps.resize;
            cone.vertices(:,2) = cone.vertices(:,2) * -1;
            cones = cone;
            for i = 1:size(sub,1)
                shiftSA = sub(i,:);
                cones(i) = cone;
                cones(i).vertices = cone.vertices + repmat(shiftSA,[size(cone.vertices,1) 1]);
                cs = cones(i);
                cs.faceColor = [.8 .8 .8]; %markerCol(m,:);
                cs.faceAlpha = 0.75;
                patch(cs)
                hold on
            end
%         end
    end
    
    %% Draw cells
    % cellList = 125;
    aL = lightangle(0,45) ;
    trackCells = [];
    
    for i = 1:length(subG)
        %subCell = names2Subs(obI,dsObj,cellList(i));
        sub = subG{i};
        if ~isempty(sub)
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
            tic
            if isempty(smallSub)
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
                %     cellDat(i).subs = sub;
                %     cellDat(i).fv = fv;
            end
        end
    end
    %aL = lightangle(0,45) ;
    
    
    %% movie
    tag = 'testMove';
    frames = 360;
    el = 0;%zeros(frames,1);
    az = 0; %0:360/frames:359;
    savefig([obMovDir tag '.fig'])
    
    imageName = sprintf('%sspringRun_%s%05.0f_01h.png',obMovDir,tag,i);
    
    startAngle = [20 20];
    
    view(startAngle)
    lightangle(aL,19.8, -22.8)
    pause(.01)
    axis off
    
    if 0
        set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
        set(gcf, 'InvertHardCopy', 'off');
        print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
    end
    
    
    %%
    
    
    for i = 1:length(az);
        view([az(i)+startAngle(1) startAngle(2)])
        lightangle(aL,az(i)+10+startAngle(1),startAngle(2)+30)
        %lightangle(aL,82.4,110)
        
        axis off
        set(gcf,'PaperUnits','points','PaperPosition',[1 1 512 512])
        %runSprings(springDat,allResults{1})
        set(gcf, 'InvertHardCopy', 'off');
        imageName = sprintf('%srot_%d_%05.0f.png',obMovDir,t,i);
        %print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
        %print(gcf, imageName, '-dpng','-opengl','-r72')
        if shouldWrite
            print(gcf,imageName,'-dpng','-r256','-opengl','-noui')
        end
        pause(.01)
    end
    disp(sprintf('finished position %d of %d',t,length(checkSyn)));
    pause(.01)
    end
end


