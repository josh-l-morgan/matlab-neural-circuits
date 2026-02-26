
clear all
%% Load data
load('MPN.mat')
if ~exist('MPN','var')
    MPN = GetMyDir;
end

synDir = [MPN 'synPos3\'];
if ~(exist(synDir,'dir')),mkdir(synDir); end

if ~exist('obI','var') | ~exist('dsObj','var')
    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
end

disp('showing')
dSamp = [8 8 4];


colMap = hsv(256);
colMap = cat(1,[0 0 0],colMap);
viewProps.dim =3;

%% Set Variables
targCells = [108];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% n popsEach possible pairing between RGCs and TCs in the subnetworks of cell A and B was classified according to whether the RGC and TC had the same seed cell association or a different seed cell association.
allEdges = obI.nameProps.edges;
preTarg = preTo(allEdges,125);
postTarg = postTo(allEdges,125);
%%Reference lists
showCell = [125 preTarg(:,1)' postTarg(:,1)'];
showCell = [125];


showCellNames = cat(2,num2cell(showCell));

col = [1 1 1];


%% Get synapse colors 

 plusOne = 125;
    %seedList = [ 108  201 109 ];
     % useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
   useList = obI2nodes_rtl_plus1(obI,plusOne);

    cellList = unique([plusOne; useList.preList(:); useList.postList(:);forceNodes]);
    
    %seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
    preSyn = allEdges(allEdges(:,1)==plusOne,:);
    [a b] = sort(preSyn(:,2),'ascend');
    preSyn = preSyn(b,:);
    
    postSyn = allEdges(allEdges(:,2) == plusOne, :);
    [a b] = sort(postSyn(:,1),'ascend');
    postSyn = postSyn(b,:);
    
    
    preNum = preTo(allEdges,plusOne);
    postNum = postTo(allEdges,plusOne);
    
    %% filter for types
    
    cellList = postNum(:,1);
    use = zeros(length(cellList),1);
    for i = 1:length(cellList)
        
        currCell = cellList(i);
        targIDX = find(obI.cell.name==currCell);
        targ  = obI.cell.mainObID(targIDX);
        
        if obI.nameProps.rgc(targ)
            use(i) = 0;
        elseif obI.nameProps.tcr(targ)
            use(i) = 1;
        elseif obI.nameProps.lin(targ)
            use(i) = 0;
        else
            use(i) = 0;
        end
                
    end
    
    postNum = postNum(use>0,:);
    [a b ]= sort(postNum(:,2),'descend');
    postNum = postNum(b,:);
    
    cumSyn = cumsum(postNum(:,2))/sum(postNum(:,2));
    plot(cumSyn)
    
    histRange = [1:20];
    histConv = histc(postNum(:,2),histRange);
    [histRange' histConv]


synCol = [0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0];
conCol = [1 2 3 4 5 6]';

Prop = postNum(:,2);
targCol = zeros(length(Prop),size(conCol,2))+.2;
targCol(Prop==1,:) = repmat(conCol(1,:),[sum(Prop==1) 1]);
targCol((Prop>1) & (Prop<=3),:) = repmat(conCol(2,:),[sum((Prop>1) & (Prop<=3)) 1]);
targCol((Prop>3) & (Prop<=6),:) = repmat(conCol(3,:),[sum((Prop>3) & (Prop<=6)) 1]);
targCol((Prop>6) & (Prop<=10),:) = repmat(conCol(4,:),[sum((Prop>6) & (Prop<=10)) 1]);
targCol(Prop>10,:) = repmat(conCol(5,:),[sum(Prop>10) 1]);


%% Get synapse positions
anchors = double(obI.colStruc.anchors);

dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;


anchors(:,1) = anchors(:,1)*dSamp(1);
anchors(:,2) = anchors(:,2)*dSamp(2);
anchors(:,3) = anchors(:,3)*dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;
synapses = obI.nameProps.edges(:,1:3);

c = 0;
useSyn = cell(max(targCol),1);
for i = 1:size(synapses,1)
    syn = synapses(i,:);
    if syn(2) == 125
        targ = find(postNum(:,1)==syn(1),1);
        if ~isempty(targ)
            c = c+1;
            g = targCol(targ);
            useSyn{g} = cat(1,useSyn{g},anchors(syn(3),:));
            
        end
    end
    
    
end







cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end





fsize = double(max(cat(1,dsObj.subs),[],1))+100;
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [0 0 0; fsize];


%% Display Variables

viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = 3;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;



degs = [0];
for d = 1:length(degs);
    
    viewProps.degRot = degs(d);
    
    
    %% Display Cells
    
    %I_topSum = showCellsAndMore(viewProps);
    I_topSum = stereoCellsAndMoreFull_PS(viewProps);
    
    %image(uint8(I_topSum*1))
    
    %% Draw synapses
    %%Ball
    ballRad = 6;
    ball = ones(ballRad*2+1,ballRad*2+1,ballRad*2+1);
    [y x z] = ind2sub(size(ball),find(ball));
    dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
    ball(:) = dists;
    
    ball1 = ball;
    ball1(ball1>ballRad) = 0;
    ballSum1 = max(ball1>0,[],3);
    
    ball2 = ball;
    ball2(ball2>5) = 0;
    ballSum2 = sum(ball2>0,3);
    ballSum2 = ballSum2 * 200/max(ballSum2(:));
    
    ballSum = ballSum1/100+ballSum2;
    
    
    dim = viewProps.dim;
    if dim == 1
        dims = [3 2];
    elseif dim == 2
        dims = [3 1];
    elseif dim == 3
        dims = [1 2];
    end
    
    %% Draw syn from track sheet
    synPos = getSynPos(1);
    
%     useSyn{1} = synPos.allPos;
%     useSyn{2} = synPos.prePos;
%     useSyn{3} = synPos.postRGCPos;
%     useSyn{4} = synPos.postUnkPos;
%     useSyn{5} = [11666 13942 3417];
%     synCol= [ 0 0 0; 2 0 0; 0 1 0; 0 0 3; 0 0 0; 0 0 0]* 1 ;
%     %synCol= [ 0 0 0; 0 1 1; 1 0 1; 0 0 0; 0 0 0; 0 0 0] ;
%     
    
    showSynTrack  = I_topSum *0;
    
    for p = 1:length(useSyn)
        anchors = double(useSyn{p});
%         anchors = anchors(:,[2 1 3]);
%         
%         dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
%         
%         anchors(:,1) = anchors(:,1)*dSamp(1);
%         anchors(:,2) = anchors(:,2)*dSamp(2);
%         anchors(:,3) = anchors(:,3)*dSamp(3);
%         
%         anchors = round(anchors);
%         anchors(anchors<1) = 1;
%         
        useAnchors = anchors;
        
        viewProps.points = useAnchors;
        useAnchors = stereoPoints_PS(viewProps);
        
        
        Isize = [viewProps.viewWindow(2,dims(1)) - viewProps.viewWindow(1,dims(1)) ....
            viewProps.viewWindow(2,dims(2)) - viewProps.viewWindow(1,dims(2))];
        
        
        anchorInd = sub2ind(Isize, useAnchors(:,dims(1)) - viewProps.viewWindow(1,dims(1)),...
            useAnchors(:,dims(2)) - viewProps.viewWindow(1,dims(1)));
        uAnchInd = unique(anchorInd);
        if length(uAnchInd)>1
            valAnch = histc(anchorInd,uAnchInd);
        else
            valAnch = 1;
        end
        
        synImage = zeros(Isize);
        synImage(uAnchInd) = valAnch;
        synImage = convn(synImage,ballSum,'same');
        
        for c = 1:3
            showSynTrack(:,:,c) =  showSynTrack(:,:,c) + synImage * synCol(p,c);
        end
        
        %image(uint8(showSynTrack*2));
        
    end
    
    maskSyn = repmat(sum(showSynTrack,3),[1 1 3]);
    
    showSynTrack(~maskSyn) = I_topSum(~maskSyn);
    %showSynTrack(~maskSyn) = I_topSum(~maskSyn) .* 0 + 255 - I_topSum(~maskSyn)*3;
    %showSynTrack = showSynTrack .* 0 + 255 - showSynTrack;
    
    image(uint8(showSynTrack*1)),pause(.1)
    
    rotDir = 'C:\Users\jlmorgan\Documents\LIN\images\rotSynDistribution5\';
    iNam = sprintf('rot%04.0f.png',d);
    if ~exist(rotDir,'dir'),mkdir(rotDir),end
    imwrite(uint8(showSynTrack ),[rotDir iNam])
    
end %end rotation

%}

%% Draw syn from anchors

anchors = double(obI.colStruc.anchors);

dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;

anchors(:,1) = anchors(:,1)*dSamp(1);
anchors(:,2) = anchors(:,2)*dSamp(2);
anchors(:,3) = anchors(:,3)*dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;
synapses = obI.nameProps.edges(:,1:3);



showSyn  = I_topSum *0;
for p = 1:length(prePop)
    
    usePre = [];
    for i = 1:length(prePop{p})
        usePre = cat(1,usePre,find(synapses(:,2)== prePop{p}(i)));
    end
    
    usePost = [];
    for i = 1:length(postPop{p})
        usePost = cat(1,usePost,find(synapses(:,1)== postPop{p}(i)));
    end
    
    foundSyn = intersect(usePre,usePost);
    
    useSynID = synapses(foundSyn,3);
    useAnchors = anchors(useSynID,:);
    
    Isize = [viewProps.fsize(dims(1)) viewProps.fsize(dims(2))];
    
    anchorInd = sub2ind(Isize, useAnchors(:,dims(1)),useAnchors(:,dims(2)));
    uAnchInd = unique(anchorInd);
    if length(uAnchInd)>1
        valAnch = histc(anchorInd,uAnchInd);
    else
        valAnch = 1;
    end
    
    %
    % for d = 1:length(anchorInd)
    %
    %
    % end
    
    synImage = zeros(Isize);
    synImage(uAnchInd) = valAnch;
    synImage = convn(synImage,ballSum,'same');
    
    for c = 1:3
        showSyn(:,:,c) =  showSyn(:,:,c) + synImage * synCol(p,c);
    end
    
    %image(uint8(showSyn*2));
    
    
end

maskSyn = repmat(sum(showSyn,3),[1 1 3]);

showSyn(~maskSyn) = I_topSum(~maskSyn) .* 0 + 255 - I_topSum(~maskSyn);
image(uint8(showSyn*1))
% rotSyn = imrotate(showSyn,pi);



% image(uint8(rotSyn)*10)

%% Write images


%{
imwrite(uint8(showSyn*.1),[cellPicDir 'syn130.png'])
imwrite(uint8(I_topSum*1),[cellPicDir 'targ108.png'])

%}

%{
cropSyn = showSyn;
[y x] = find(sum(cropSyn,3)>0);
cropSyn = cropSyn(min(y):max(y),min(x):max(x),:);
image(uint8(cropSyn * 4)),pause(.01)
%}

%{

imwrite(uint8(cropSyn),sprintf('%ssynPos_%06.0f.png',synDir,showCell))


imwrite(uint8(cropSyn),sprintf('%ssynPosNoSyn_%06.0f.png',synDir,showCell))

%}

