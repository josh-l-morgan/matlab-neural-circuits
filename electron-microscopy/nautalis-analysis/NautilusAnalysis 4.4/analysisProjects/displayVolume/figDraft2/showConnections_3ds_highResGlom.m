
clear all
%% Load data
%load('MPN.mat')
MPN = 'D:\LGNs1\Segmentation\VAST\S4\joshm\export_target\glomA\exportGlomA_14+08+25_fullREs_mat\'


figDir = [MPN 'figDraft1\'];
if ~(exist(figDir,'dir')),mkdir(figDir); end

if ~exist('obI','var') | ~exist('dsObj','var')
    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
end

disp('showing')
dSamp = [2 2 1];

colMap = hsv(256);

viewProps.dim =3;

%% Set Variables
targCells = [ 108 201  903 907 ];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% n pops

if 1
   
    useCells = [conTo(2).rgcList];
    %useCells = useCells(1:3);
    
    col = colMap(ceil((1:length(useCells))*256/length(useCells)),:);
    col = col(randperm(size(col,1)),:);
    targCol = [1 .5 1; .5 1 1];
    col = [targCol; col];
    
    
     load('.\data\colormaps\col_glomA1.mat')
    % save('.\data\colormaps\col_glomA1.mat','col')
    col(1:2,:) = [1 .7 1; .7 1 1] * .8;
    
    showCellNames = cat(2,num2cell([ 201 203 useCells]));

end


%% prep image
cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = double([1 1 1; fsize(1) fsize(2) fsize(3)]);
viewProps.viewWindow = double([1 1 3000; fsize(1) fsize(2) fsize(3)]);

%% Display Variables

viewProps.maxScaleFactor = 1;
viewProps.sumScaleFactor = 0;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;
viewProps.resamp = [.75 .5 1];



%% Display Cells
if 0
    I_topSum = showCellsAndMore(viewProps);
    image(uint8(I_topSum*1))
end

if 0
    Ist = showStereo(viewProps);
    image(uint8(Ist*.1) )
    % 
end


if 1 %% render objects
    objDir = 'D:\LGNs1\Analysis\movies\glomA\obj7\'
    if ~exist(objDir,'dir'),mkdir(objDir), end
    showSTL(viewProps,objDir)
end


if 0
   %% 
    viewProps.keepRat = 1;
    viewProps.contrast = 1;
    viewProps.gamma = .9;
    viewProps.dilate = 0;
    
    rotDir = 'D:\LGNs1\Analysis\movies\glomA\rot5\'
    if ~exist(rotDir,'dir'),mkdir(rotDir), end
    
    degs = [0:1:359];
    %degs = [-1.5 1.5];
    %degs = [0:1:359];
    %degs = [0:3:15];
    
    I_stack = showRotationPerspective(viewProps,degs,rotDir);
    %playRot(I_stack);
   %% 
end

if 0
    %%
       % I = Ist;
        I = I_topSum;
    viewProps.keepRat = 1;
    viewProps.contrast = 1;
    viewProps.gamma = .9;
    viewProps.dilate = 0;

    
        I = tweakI(I,viewProps);
        
%             Iblank = sum(I,3) == 0;
% 
%          I(repmat(Iblank,[1 1 3]) > 1) = 255;

    image(uint8(I))
        
    %%
    %   imwrite(uint8(I),[cellPicDir 'ABmix.png'])
    %%  imwrite(uint8(I),[objDir 'exampleGlomA.png'])
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

