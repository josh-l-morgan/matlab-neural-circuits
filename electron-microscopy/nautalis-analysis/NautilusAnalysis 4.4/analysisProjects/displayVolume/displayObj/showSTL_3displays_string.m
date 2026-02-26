
clear all
%% Load data
load('MPN.mat')
if ~exist('MPN','var')
    MPN = GetMyDir
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
%colMap = cat(1,[0 0 0],colMap);

viewProps.dim =3;

%% Set Variables
targCells = [ 108 201  903 907 ];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% n pops

if 1
    
    useCells = cat(1,conTo([1:4]).tcrList);
    showCellNames = cat(2,num2cell(useCells));
    col = getList_seedColor(targCells,useCells);
        
end

if 0
    
    if 0
        
        load('.\data\clade\cladePick_six2.mat')
        
        useCells = [conTo([1:4]).tcrList];
        
        
        colorTable = [ 1 0 0; 0 1 0; 0 1 1.5; 1 0 1.5]
        %colorTable = [ 0 .2 1; 1 .2 0; 0 1 1.5; 1 0 1.5]*2
        col = [];
        allMembers = [];
        useClade = [ 1 2 3 4];
        
        for i = 1:length(useClade)
            members = conTo(useClade(i)).targ;
            %members = getList_glomID(useClade(i));
            %members = intersect(useCells,members);
            members = members(:)';
            allMembers = [allMembers members];
            col = cat(1,col,repmat(colorTable(useClade(i),:),[length(members),1]));
            
        end
        
        %%Reference lists
        showCell = [allMembers ];
        col1 = cat(1,col);
        showCellNames1 = cat(2,num2cell(showCell));
        
    else
        
        col1 = [];
        showCellNames1 = {};
    end
    
    if 0
        
        load('.\data\clade\cladePick_six2.mat')
        
        useCells = [conTo([1:4]).tcrList];
        
        
        colorTable = [ 1 0 0; 0 1 0; 0 1 1.5; 1 0 1.5]
        %colorTable = [ 0 .2 1; 1 .2 0; 0 1 1.5; 1 0 1.5]*2
        col = [];
        allMembers = [];
        useClade = [ 1 2 3 4];
        
        for i = 1:length(useClade)
            members = conTo(useClade(i)).tcrList;
            %members = getList_glomID(useClade(i));
            %members = intersect(useCells,members);
            members = members(:)';
            allMembers = [allMembers members];
            col = cat(1,col,repmat(colorTable(useClade(i),:),[length(members),1]));
            
        end
        
        
        %%Reference lists
        showCell = [allMembers ];
        col1 = cat(1,col);
        showCellNames1 = cat(2,num2cell(showCell));
        
    else
        %
        %     col1 = [];
        %     showCellNames1 = {};
    end
    
    
    if 0
        
        checkGlom = getList_glomID(1);
        showCell = [108 201 903 907];
        showCellNames2 = cat(2,num2cell(showCell),'axon 108','axon 201','axon 903', 'axon 907');
        col2 = [ 1 0 0; 0 1 0; 0 1 1.5; 1 0 1.5; 2 2 0 ;2 2 0 ; 2 2 0 ;2 2 0 ] * .2 + .9;
        
    else
        
        col2 = [];
        showCellNames2 = {};
        
    end
    
    showCellNames = cat(2,showCellNames1,showCellNames2);
    col = cat(1,col1,col2)
end





%% prep image
cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end



fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];
%}


%% Select Window

%{

mid = [11422, 19706, 4689]; % glomB
mid = [14815, 18606, 2581] %glomG
mid = [12492, 21454, 5875] %massive input on 30001 225
mid = [13496, 19858, 3743] %glomA
mid = [11375, 15869, 4001] %giant diverge 9101 onto 903 and 919

mid = round([mid(1)/8 mid(2)/8 mid(3)/4]);
rad = [50 50 50];
viewProps.viewWindow = [mid(2) - rad(1)  mid(1) - rad(2)  mid(3) - rad(3) ; ...
    mid(2) + rad(1)  mid(1) + rad(2)  mid(3) + rad(3) ];



% viewProps.viewWindow = [2300  1400  600 ; ...
%     2400  1600  1200 ];

%Target 2350 1500 700, Y 18800, x 1200, z = 2800

%}
viewProps.viewWindow = double([1 1 1; fsize]);
viewProps.viewWindow = double([fsize(1)/3 1 1; fsize(1) fsize(2) fsize(3)]);

%% Display Variables

viewProps.maxScaleFactor = 5;
viewProps.sumScaleFactor = 0;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;



%% Display Cells
if 1
    I_topSum = showCellsAndMore(viewProps);
    image(uint8(I_topSum*1))
end

if 0
    Ist = showStereo(viewProps);
    image(uint8(Ist) )
    
end


if 0
    rotDir = 'D:\LGNs1\Analysis\movies\fourSeeds\seed4fullTCR_saturate\'
    if ~exist(rotDir,'dir'),mkdir(rotDir), end
    degs = [90];
    degs = [0:1:359];
    %degs = [0:3:15];
    I_stack = showRotationPerspective(viewProps,degs,rotDir);
  
    %playRot(I_stack);
    
end

if 0
    %%
    I = I_topSum;
    gamma = .2;
    I = I.^gamma * 256/256^gamma;
    SE = strel('disk',1);
    %I = imdilate(I,SE);
    I = I * 1;
    keepRat = .2;
    %I = uint8keepRat(I*.5)*keepRat + uint8(I*.5)*(1-keepRat);
    image(uint8(I))
    %%
    %   imwrite(uint8(I),[cellPicDir 'fourSeedTCR_gamma1.png'])
    
end

if 0
    
    %% Pick synapses
    cellList = obI.cell.name;
    clear prePop postPop popCol
    for i = 1:length(useClade)
        prePop{i} = conTo(useClade(i)).rgcList;
        postPop{i} = conTo(useClade(i)).tcrList;
        popCol{i} = colorTable(useClade(i),:);
        
    end
    
    
    
    %% Draw synapses
    %%Ball
    ballRad = 6;
    ball = ones(ballRad*2+1,ballRad*2+1,ballRad*2+1);
    [y x z] = ind2sub(size(ball),find(ball));
    dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
    ball(:) = dists;
    ball(ball>ballRad) = 0;
    ballSum = sum(ball>0,3);
    ballSum = ballSum * 200/max(ballSum(:));
    
    anchors = double(obI.colStruc.anchors);
    
    dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
    
    anchors(:,1) = anchors(:,1)*dSamp(1);
    anchors(:,2) = anchors(:,2)*dSamp(2);
    anchors(:,3) = anchors(:,3)*dSamp(3);
    anchors = round(anchors);
    anchors(anchors<1) = 1;
    synapses = obI.nameProps.edges(:,1:3);
    
    
    dim = viewProps.dim;
    if dim == 1
        dims = [3 2];
    elseif dim == 2
        dims = [3 1];
    elseif dim == 3
        dims = [1 2];
    end
    
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
            showSyn(:,:,c) =  showSyn(:,:,c) +synImage*popCol{p}(c) ;
        end
        image(uint8(showSyn*2));
        
        
    end
    
    
    maskSyn = repmat(sum(showSyn,3),[1 1 3]);
    showBoth = showSyn;
    showBoth(~maskSyn) = I_topSum(~maskSyn)/2;
    SE = strel('disk',1);
    Idia = imdilate(showSyn,SE);
    %Idia = imdilate(Idia,SE);
    image(uint8(Idia*.5))
    % rotSyn = imrotate(showSyn,pi);
    % image(uint8(rotSyn)*10)
    
end
%% Write images


%{
imwrite(uint8(showSyn),[cellPicDir 'CBsyn.png'])
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

