
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





%% Get synapse colors


tri = getTriads;

prim125 = find(tri.triCell(:,1)==125);
useTri = prim125;
for u = 1:length(prim125);%1:size(tri.triCell,1);
    useTri = prim125(u);
    
    triPoints = drawTriads(tri,useTri);
    
    synCol = cat(1,triPoints.lineCol,triPoints.ballCol);
    synRad = cat(1,triPoints.lineRad,triPoints.ballRad);
    synType = cat(1,triPoints.lineRad*0+1,triPoints.ballRad*0+2);
    useSyn = cat(2,triPoints.lineGroup, triPoints.ballGroup);
    primCell = unique(triPoints.cellGroup(:,1));
    secCell = unique(triPoints.cellGroup(:,2));
    tertCell = unique(triPoints.cellGroup(:,3));
    
    
    
    showCellNames = cat(2,num2cell([primCell; secCell; tertCell]));
    
    col = [repmat([1 0 0],[length(primCell) 1]);...
        repmat([0 1 0],[length(secCell) 1]); repmat([0 0 1],[length(tertCell) 1])];
    
    
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
    
    
    
    degs = [0:90 ];
    for d = 1:length(degs);
        
        viewProps.degRot = degs(d);
        
        
        %% Display Cells
        
        %I_topSum = showCellsAndMore(viewProps);
        I_topSum = stereoCellsAndMoreFull_PS(viewProps);
        
        %image(uint8(I_topSum*1))
        
        %% Draw synapses
        
        dim = viewProps.dim;
        if dim == 1
            dims = [3 2];
        elseif dim == 2
            dims = [3 1];
        elseif dim == 3
            dims = [1 2];
        end
        
        %% Draw syn from track sheet
        %synPos = getSynPos(1);
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
            useAnchors = anchors;
            
            if synType(p) == 1
                ballSum = ballPic(synRad(p));
            else
                ballSum = ringPic(synRad(p));
            end
            
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
        
        
        rotDir = 'C:\Users\jlmorgan\Documents\LIN\images\rotSynDistribution17\';
        %     rotFold = sprintf('%striad_%04.0f\\',rotDir,useTri);
        iNam = sprintf('tri%04.0f_rot%04.0f.png',useTri,d);
        if ~exist(rotDir,'dir'),mkdir(rotDir),end
        writeImage = uint8(showSynTrack(1000:2500,500:2000,:));
        imwrite(uint8(writeImage ),[rotDir iNam])
        
    end %end rotation
    
end %end tri
%{


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
%}


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
