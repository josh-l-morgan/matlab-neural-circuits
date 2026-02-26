clear all
load('MPN.mat');
load([MPN 'obI.mat']);
targCell = 125;
maxDist = 5;
iptsetpref('ImshowBorder','tight');

disp('combining spreadsheet and segmentation synapses')
sm = addDatToSynMat(obI)
disp('finding topological distance between synapses on target cell within max range')
sm  = getTopoEucDistBetweenSyn(sm);
sm = getTopoEucDistBetweenSkelAndSyn(sm)

%% analyze Motifs


isTarg =(sm.pre == targCell) | (sm.post == targCell);
isRGC = sm.preClass == 1;
isTC = sm.postClass == 2;
isPreLIN = sm.preClass == 3;
isPostLIN = sm.postClass == 3;
isPreUNK = sm.preClass == 4;


%search all RGC input, check for non RGC TC
%search all LIN input, check for isolated LIN output
%search for all UNK input


%% find checkIn relative motifs



[checkRgc] = checkMotifAroundSyns(sm,find(isRGC & isTarg),maxDist);
[checkPreLin] = checkMotifAroundSyns(sm,find(isPreLIN & (sm.post == targCell)),maxDist);
[checkUnk] = checkMotifAroundSyns(sm,find(isPreUNK & isTarg),maxDist);
[checkPostLin] = checkMotifAroundSyns(sm,find(isPostLIN & (sm.pre == targCell)),maxDist);
[checkTc] = checkMotifAroundSyns(sm,find(isTC & isTarg),maxDist);



motifCodes = [...
    '1 = RGC only, 2 = TC only, 3 = LINin only, 4 = LINout only, 5 = UNKin only, ' ...
    '6 = recip, 7 = autapse,'...
    '8 = RGC-LIN-TC diad, 9 = RGC-LIN-TC triad, 10 = RGC-LIN-TC mult triad, ' ...
    '11 = RGC-LIN-LIN diad, 12 = RGC-LIN-LIN IN triad,  13 = RGC-LIN-LIN out triad, 14 = RGC-LIN-LIN recip triad, ' ...
    '15 = LIN-LIN-TC diad (top), 16 = LIN-LIN-TC diad (mid), '...
    '17 = LIN-LIN-TC triad(top), 18 = LIN-LIN-TC triad(mid), 19 = LIN-LIN-TC recip triad, '...
    '20 = LLL diad, 21 = LLL tri top, 22 = LLL tri mid, 23 = LLL tri bot, 24 = LLL recip triad, '...
    '25 = UNK-LIN-TC(diad), 26 = UNK-LIN-TC triad'];



%% classify motifs


for i = 1:30, mPos{i} = []; end

for i = 1:length(checkRgc.m)
    m = checkRgc.m(i);
    if m.preOnly
        mPos{1} = [mPos{1};checkRgc.pos(i,:)];
    elseif m.preTcTri
        mPos{9} = [mPos{9};checkRgc.pos(i,:)];
    elseif m.preTcDi
        mPos{8} = [mPos{8};checkRgc.pos(i,:)];
    elseif m.preLinDi
        mPos{11} = [mPos{11};checkRgc.pos(i,:)];
    elseif m.preLinInTri
        mPos{13} = [mPos{13};checkRgc.pos(i,:)];
    elseif m.preLinOutTri
        mPos{12} = [mPos{12};checkRgc.pos(i,:)];
    end
end

for i = 1:length(checkPreLin.m)
    m = checkPreLin.m(i);
    if m.preOnly
        mPos{3} = [mPos{3};checkPreLin.pos(i,:)];
    elseif m.preTcTri
        mPos{18} = [mPos{18};checkPreLin.pos(i,:)];
    elseif m.preTcDi
        mPos{16} = [mPos{16};checkPreLin.pos(i,:)];
    elseif m.preLinDi
        mPos{20} = [mPos{20};checkPreLin.pos(i,:)];
    elseif m.preLinInTri
        mPos{23} = [mPos{23};checkPreLin.pos(i,:)];
    elseif m.preLinOutTri
        mPos{22} = [mPos{23};checkPreLin.pos(i,:)];
    end
end

for i = 1:length(checkUnk.m)
    m = checkUnk.m(i);
    if m.preOnly
        mPos{5} = [mPos{5};checkUnk.pos(i,:)];
    elseif m.preTcTri
        mPos{26} = [mPos{26};checkUnk.pos(i,:)];
    elseif m.preTcDi
        mPos{25} = [mPos{25};checkUnk.pos(i,:)];
    elseif m.preLinDi
        mPos{27} = [mPos{27};checkUnk.pos(i,:)];
    elseif m.preLinInTri
        mPos{28} = [mPos{28};checkUnk.pos(i,:)];
    elseif m.preLinOutTri
        mPos{29} = [mPos{29};checkUnk.pos(i,:)];
    end
end

for i = 1:length(checkPostLin.m)
    m = checkPostLin.m(i)
    
    if m.preOnly
        mPos{4} = [mPos{4};checkPostLin.pos(i,:)];
    elseif m.preLinDi
        mpfos{15} = [mPos{15};checkPostLin.pos(i,:)];
    elseif m.preTcTri
        mPos{17} = [mPos{17};checkPostLin.pos(i,:)];
        
    elseif m.recTcTri
        mPos{19} = [mPos{19};checkPostLin.pos(i,:)];
    elseif m.recRgcTri
        mPos{14} = [mPos{14};checkPostLin.pos(i,:)];
    elseif m.recLinInTri
        mPos{24} = [mPos{24};checkPostLin.pos(i,:)];
    elseif m.recLinOutTri
        mPos{24} = [mPos{24};checkPostLin.pos(i,:)];
    elseif m.recLinInOut
        mPos{24} = [mPos{24};checkPostLin.pos(i,:)];
    elseif m.recLinOutIn
        mPos{24} = [mPos{24};checkPostLin.pos(i,:)];
    end
    if m.reciprocal
        mPos{6} = [mPos{6};checkPostLin.pos(i,:)];
    end
    if m.autapse
        mPos{7} = [mPos{7};checkPostLin.pos(i,:)];
    end
end

for i = 1:length(checkTc.m)
    m = checkTc.m(i);
    if m.postOnly
        mPos{2} = [mPos{2};checkTc.pos(i,:)];
    end
end


%%find clusters

%distances between skeleton points
sPos = sm.skelPos;
sDist = sqrt((sPos(:,1)-sPos(:,1)').^2 +  (sPos(:,2)-sPos(:,2)').^2 + ...
    (sPos(:,3)-sPos(:,3)').^2 ); % distance to other points on arbor
minClust = 3;
clustDist = 5;
clear clust
for i =  1:length(mPos)
    pos = mPos{i};
    
    if ~isempty(pos)
        dist = sqrt((pos(:,1)-sPos(:,1)').^2 +  (pos(:,2)-sPos(:,2)').^2 + ...
            (pos(:,3)-sPos(:,3)').^2 ); %distance to motif
        nearNum = sum(dist<=clustDist,1);
        eraseNum = nearNum;
        keepNum = nearNum * 0;
        c = 0;
        while sum(eraseNum)
            c = c + 1;
            targMax = find(eraseNum == max(eraseNum),1);
            keepNum(targMax) = nearNum(targMax); %record peak
            eraseNum(sDist(:,targMax)<(clustDist*2)) = 0; %erase nearby
        end
        isClust = find(keepNum>=minClust);
        clust(i).pos = sPos(isClust,:);
        clust(i).val = keepNum(isClust);
        
    else
        clust(i).pos = [];
        clust(i).val = [];
    end
    
end


%%find clustered

%distances between skeleton points
sPos = sm.skelPos;
sDist = sqrt((sPos(:,1)-sPos(:,1)').^2 +  (sPos(:,2)-sPos(:,2)').^2 + ...
    (sPos(:,3)-sPos(:,3)').^2 ); % distance to other points on arbor
minClust = 3;
clustDist = 5;
for i =  1:length(mPos)
    pos = mPos{i};
    if ~isempty(pos)
        dist = sqrt((pos(:,1)-pos(:,1)').^2 +  (pos(:,2)-pos(:,2)').^2 + ...
            (pos(:,3)-pos(:,3)').^2 );
        %distance to motif
        nearNum = sum(dist<=clustDist,1);
        mClust{i} = nearNum>=(minClust );
    else
        mClust{i} = [];
    end
end

%% Show result


% save(['D:\LGNs1\Analysis\motifImages\listMotifWorkSpace2.mat'])
% load(['D:\LGNs1\Analysis\motifImages\listMotifWorkSpace2.mat'])

oneAtATime = 1;
shouldWrite = 1;
motNum = 30;
skelWidth = 4;
circSize = 80^2;
circWidth = 4;
skelCol = [1 0 0];
markWidth  = 5;

obMovDir = 'D:\LGNs1\Analysis\motifImages\batch12rev\'
if ~exist(obMovDir,'dir'), mkdir(obMovDir);end

motText = {'1r' '1t' '1li' '1lo' '1u' '2b' '2a' '3d' '3t' '3m'...
    '4d' '4t' '4t'  '4r'  '5d' '5d' '5t' '5t' '5r' '5d'...
    '6t' '6m' '6b'  '6r'  '7d' '7t' '7l' '7l' '7l' '0'};

if oneAtATime == 1
    motCol = {'w' 'w' 'w' 'w' 'w' 'w' 'w' 'w' 'w' 'w'...
        'w' 'w' 'w'  'w'  'w' 'w' 'w' 'w' 'w' 'w'...
        'w' 'w' 'w'  'w'  'w' 'w' 'w' 'w' 'w' 'w'};
    for i = 1:motNum, motCol{i} = 'r'; end
    for i = 1:motNum, clustCol{i} = [0 0 1]; end
    
    motMark = {'p' 'p' 'p' 'p' 'p' 'p' 'p' 'p' 'p' 'p'...
        'p' 'p' 'p'  'p'  'p' 'p' 'p' 'p' 'p' 'p'...
        'p' 'p' 'p'  'p'  'p' 'p' 'p' 'p' 'p' 'p'};
    for i = 1:motNum, motMark{i} = 'o'; end
    
    motSize = [200 200 200 200 200 200 200 200 200 200 ...
        200 200 200  200  200 200 200 200 200 200 ...
        200 200 200  200  200 200 200 200 200 200];
    for i = 1:motNum, motSize(i) =40^2; end
    
    for i = 1:motNum, motSize(i) = 30^2; end
    motSize(7) = 40^2;
    for i = 1:motNum, motCol{i} = [.5 1 1]; end
    for i = 1:motNum, clustCol{i} = [0 0 1]; end
    for i = 1:motNum, edgeCol{i} = [0 0 0]; end
    
    showClust = [];
    
elseif oneAtATime == 0
    
    motCol = {'g' 'b' 'r' 'r' [1 .5 0] 'r' 'w' 'c' 'c' 'c'...
        'y' 'y' 'y'  'y'  'm' 'm' 'm' 'm' 'm' 'm'...
        'r' 'r' 'r'  'r'  [1 .5 1] [1 .5 1] [1 .5 0] [1 .5 0] 'm' 'w'};
    clustCol = motCol;
    
    motMark = {'o' 'o' 'o' 'o' 'o' 'p' 'p' 's' '^' '^'...
        's' '^' '^'  '^'  's' 's' 'v' '^' '^' 's'...
        'v' '^' '^'  '^'  's' '^' '^' '^' '^' 'o'};
    
    
    motSize = [50 50 50 50 50 90 1000 60 100 100 ...
        60 100 100  100  60 60 200 100 100 60 ...
        200 100 100  100  60 100 100 100 100 10];
    
    showClust = [];
    
elseif oneAtATime == 2;
    
    motCol = {'g' 'b' 'r' 'r' [1 .5 0] 'r' 'w' 'c' 'c' 'c'...
        'y' 'y' 'y'  'y'  'm' 'm' 'm' 'm' 'm' 'm'...
        'r' 'r' 'r'  'r'  [1 .5 1] [1 .5 1] [1 .5 0] [1 .5 0] 'm' 'w'};
    clustCol = motCol;
    
    motMark = {'o' 'o' 'o' 'o' 'o' 's' 'p' 's' '^' '^'...
        's' '^' '^'  '^'  's' 's' '^' '^' '^' 's'...
        '^' '^' '^'  '^'  's' '^' '^' '^' '^' 'o'};
    
    for i = 1:motNum, motSize(i) = 30^2; end
    motSize(7) = 70^2;
    for i = 1:motNum, motCol{i} = [.5 1 1]; end
    for i = 1:motNum, clustCol{i} = [0 0 1]; end
    for i = 1:motNum, edgeCol{i} = [0 0 0]; end
    
    showClust = [];
    
    motGroup = [8 9]; %RGC
     motGroup = [25 26]; % UNK
     motGroup = [16 18]; %Lin Lin
     motGroup = [15 17]; %Lin Lin Targ on top
     motGroup = [6 7]; %Lin Lin Targ on top
     motGroup = [11 12]; %rgc Lin Targ on top
    
end



hold off
% scatter3(sm.skelNodes(:,1),sm.skelNodes(:,2),sm.skelNodes(:,3),10,'filled',...
%     'markerfacecolor',[0 0 0]);

sPos = sm.skelPos;
edges = sm.skelEdges;
% plot3([sPos(edges(:,1),1) sPos(edges(:,2),1)]',...
%     [sPos(edges(:,1),2) sPos(edges(:,2),2)]',...
%     [sPos(edges(:,1),3) sPos(edges(:,2),3)]','color',skelCol,...
%     'linewidth',skelWidth);

cellSurf = 'D:\LGNs1\mergeSeg_mat\FVFiles\Patch_dSamp4_cell_125.mat';
cellSurf = 'D:\LGNs1\mergeSeg_mat\FVFiles\Patch_dSamp2_cell2_125.mat';

load(cellSurf)
fv.vertices = fv.vertices(:,[2 1 3])/2.5;
hold off
scatter3([],[],[])
patch(fv)


axis vis3d;
daspect([1,1,1])
axis tight
clight = camlight
lighting gouraud
set(gca,'Ydir','reverse')
hold on
showMot = [1]

if oneAtATime == 2
    runPos = motGroup;
else
    runPos = [1:length(mPos)];
end

for r =  1:length(runPos)
    disp(sprintf('showing %d of %d',r,length(runPos)))
    i = runPos(r);
    if oneAtATime == 1
        hold  off
        scatter3([],[],[])
        
        %         plot3([sPos(edges(:,1),1) sPos(edges(:,2),1)]',...
        %             [sPos(edges(:,1),2) sPos(edges(:,2),2)]',...
        %             [sPos(edges(:,1),3) sPos(edges(:,2),3)]','color',skelCol,...
        %             'linewidth',skelWidth);
        patch(fv)
        axis vis3d;
        daspect([1,1,1])
        axis tight
        clight = camlight
        lighting gouraud
        
        hold on
        set(gca,'Ydir','reverse')
        
    end
    
    tag = motText{i};
    pos = mPos{i};
    if ~isempty(pos)
        if sum([6 7] == i)
            use = (mClust{i} * 0)>0;
        else
            use = mClust{i};
        end
        scatter3(pos(use,1),pos(use,2),pos(use,3),motSize(i),'marker',motMark{i},...
            'MarkerFaceAlpha',.8 ,'MarkerEdgeAlpha',.8,...
            'markerfacecolor',motCol{i}, 'markeredgecolor',edgeCol{i},'linewidth',markWidth)
        scatter3(pos(~use,1),pos(~use,2),pos(~use,3),motSize(i),'marker',motMark{i},...
            'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.8,...
            'markerfacecolor',[.5 .5 .5],'markeredgecolor','k','linewidth',markWidth)
        %         scatter3(pos(~use,1),pos(~use,2),pos(~use,3),motSize(i),'marker',motMark{i},...
        %             'markerfacecolor',[.5 .5 .5],'markeredgecolor','k','linewidth',markWidth)
    end
    
    pos = clust(i).pos;
    if ~isempty(pos) & sum(showClust==i)
        scatter3(pos(:,1),pos(:,2),pos(:,3),circSize,...
            'markeredgecolor',clustCol{i},'linewidth',circWidth)
    end
    
    
    %%Show figure
    %[an el ] = view(gca);
    
    an = -115.8;
    el = -15.6;
    lightangle(clight,122,-44)
    set(gca,'view',[an el])
    set(gcf,'color',[1 1 1])
    set(gca,'color',[1 1 1])
    
    axis off
    ax = gca;
    ax.Clipping = 'off';
    set(gcf,'units','pixels','outerposition',[0 0 1000 1000]);
    %set(gca,'box','off')
    set(gca,'position',[-.00 -.05 1.15 1.15])
    
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 1024 1024])
    set(gcf, 'InvertHardCopy', 'off');
    
    pause(.01)
    if (oneAtATime == 1) & shouldWrite
        text(640,20,num2str(i),'color',[.8 .8 .8])
        title(num2str(i),'color','k')
        imageName = sprintf('%smot%d.png',obMovDir,i);
        print(gcf,imageName,'-dpng','-r512','-opengl','-noui')
    end
    
    
end

if (oneAtATime == 0) & shouldWrite
    imageName = sprintf('%sallMot.png',obMovDir);
    print(gcf,imageName,'-dpng','-r128','-opengl','-noui')
elseif (oneAtATime == 2) & shouldWrite
    imageName = sprintf('%smot_%d_%d.png',obMovDir,motGroup(1),motGroup(2));
    print(gcf,imageName,'-dpng','-r128','-opengl','-noui')
end





%
%
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = ti(1) ;
% bottom = ti(2) ;
% ax_width = ti(3) ;
% ax_height = ti(4);
% ax.Position = [left bottom ax_width ax_height];






