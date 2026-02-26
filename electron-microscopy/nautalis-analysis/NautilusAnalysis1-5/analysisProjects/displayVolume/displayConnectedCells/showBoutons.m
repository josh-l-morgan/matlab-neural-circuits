

%% Load data
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

allEdges = obI.nameProps.edges(:,[2 1]);

colMap = hsv(256);
colMap = cat(1,[0 0 0],colMap);

clear viewProps
%% Set Variables
targCells = [201];
useList = obI2cellList_seedInput(obI,targCells);



%% Find connectivity
conTo = makeConTo(obI,targCells);

%% pick cells

%%Reference lists

axnet1 =  [1000:1500];
axnet2 = [2000:2031 2038:2300];
axnet3 = [2032:2037];
oneAx = 2040;

showCell = [targCells];

%showCell = [targCells conTo(1).rgcList];

showCellNames = cat(2,num2cell(showCell));

viewProps.dim =1;
col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);
col(1,:) = 1;


%% two pops

%{
pop1 = [conTo(1).tcrList ];
pop2 = [conTo(2).tcrList ];
pop3 = [targCells];
showCellNames = num2cell([pop1 pop2 pop3]);
col = zeros(length(pop1)+length(pop2) + length(pop3),3);
col(1:length(pop1),1) = 1;
col(length(pop1)+ 1 : length(pop1) + length(pop2),2) = 1;
col(length(pop1) + length(pop2)+ 1:end,3) = 1;


%col = colMap(ceil((1:length(dsObj))*256/length(dsObj)),:);

%}


%{
%% labeled CB ?
pop1 = [ 108 201 109];
pop2 = cat(2,conTo(1).tcrList, conTo(2).tcrList);
cbList = [pop1 pop2];
cbRad =  [ones(1,length(pop1)) * 8 ones(1,length(pop2)) * 5];
cbCol =   [ones(1,length(pop1)) * 8 ones(1,length(pop2)) * 5];
%}


%% Pick synapses
prePop{1} = [preTo(allEdges,targCells)];
postPop{1} =  [targCells];

% 
% prePop{3} = [];
% postPop{3} =  [];
% 
% prePop{1} = [];
% postPop{1} =  [];


cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

cellList = obI.cell.name;




fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];
%}


%% two pops
%{
pop1 = [conTo(1).tcrList];
pop2 = [conTo(2).tcrList ];
pop3 = [conTo(3).tcrList ];
showCellNames = num2cell([pop1 pop2 pop3]);
col = zeros(length(pop1)+length(pop2) + length(pop3),3);
col(1:length(pop1),1) = 1;
col(length(pop1)+ 1 : length(pop1) + length(pop2),2) = 1;
col(length(pop1) + length(pop2)+ 1:end,3) = 1;

%}

%col = colMap(ceil((1:length(dsObj))*256/length(dsObj)),:);
% 

%% Select Window

%{

mid = [11422, 19706, 4689]; % glomB
mid = [14815, 18606, 2581] %glomG
mid = [13496, 19858, 3743] %glomA
mid = [12492, 21454, 5875] %massive input on 30001 225

mid = round([mid(1)/8 mid(2)/8 mid(3)/4]);
rad = [300 300 300];
viewProps.viewWindow = [mid(2) - rad(1)  mid(1) - rad(2)  mid(3) - rad(3) ; ...
    mid(2) + rad(1)  mid(1) + rad(2)  mid(3) + rad(3) ];
%}


% viewProps.viewWindow = [2300  1400  600 ; ...
%     2400  1600  1200 ];

%Target 2350 1500 700, Y 18800, x 1200, z = 2800

%% Display Variables

viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = 3;
viewProps.viewWindow = double([1 1 1; fsize]);
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;



%% Display Cells

I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum*1))


%% Draw boutons


viewProps.cellId = num2cell(useList.preList);
viewProps.cellId = num2cell(prePop{1});


col2 = colMap(ceil((1:length(viewProps.cellId))*256/length(viewProps.cellId)),:);
col2 = col2(randperm(size(col2,1)),:);
%col2 = [1 0 0; 0 1 1; 0 1 1; 0 1 1];
viewProps.col = col2;

I_axons = showCellsAndMore(viewProps);
image(uint8(I_axons))

%% Ball
ballRad = 32;
ball = ones(ballRad*2+1,ballRad*2+1,ballRad*2+1);
[y x z] = ind2sub(size(ball),find(ball));
dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
ball(:) = dists;
ball(ball>ballRad) = 0;
ballSum = sum(ball>0,3);
ballSum = ballSum * 200/max(ballSum(:));
[y x z ] =ind2sub(size(ball),find(ball>0));
ballSub = [y x z] - ballRad - 1;


anchors = double(obI.colStruc.anchors);
anchors(:,1) = anchors(:,1)/dSamp(1);
anchors(:,2) = anchors(:,2)/dSamp(2);
anchors(:,3) = anchors(:,3)/dSamp(3);
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

    usePre = [];
    for i = 1:length(prePop{1})
        usePre = cat(1,usePre,find(synapses(:,2)== prePop{1}(i)));
    end
        
    usePost = [];
    for i = 1:length(postPop{1})
        usePost = cat(1,usePost,find(synapses(:,1)== postPop{1}(i)));
    end
        
    foundSyn = intersect(usePre,usePost);
    
    useSynID = synapses(foundSyn,3);
    useAnchors = anchors(useSynID,:);
    
    dilateAnchors = zeros(size(useAnchors,1)*size(ballSub,1),3);
    S = size(ballSub,1);
    for a = 1:size(useAnchors,1)
        L = (a-1) * S;
        dilateAnchors(L+1:L+S,1) = ballSub(:,1) + useAnchors(a,1);
        dilateAnchors(L+1:L+S,2) = ballSub(:,2) + useAnchors(a,2);
        dilateAnchors(L+1:L+S,3) = ballSub(:,3) + useAnchors(a,3);
    end
    
    %%Get unique
    dilateAnchors(dilateAnchors<1) = 1;
    maxD = max(dilateAnchors,[],1);
    DAind = sub2ind(maxD,dilateAnchors(:,1),dilateAnchors(:,2),dilateAnchors(:,3));
    uDAind = unique(DAind);
    [y x z] = ind2sub(maxD,uDAind);
    dilateAnchors = [y x z];
    
    viewProps.mask = dilateAnchors;
    
    
    
I_boutons = showCellsAndMore(viewProps);
image(uint8(I_boutons*3))

%%
I_pre = I_boutons;
 maskSyn = repmat(sum(I_pre,3),[1 1 3]);

showSyn = I_pre;
showSyn(~maskSyn) = (I_topSum(~maskSyn))*.4;
%showSyn = I_topSum*.1;

image(uint8(showSyn))
% rotSyn = imrotate(showSyn,pi);
% image(uint8(rotSyn))

%% Write images


if 0
    %%
    
imwrite(uint8(showSyn),sprintf('%sboutons%d.png',cellPicDir,targCells))
imwrite(uint8(I_axons),sprintf('%saxons%d.png',cellPicDir,targCells))

end
%}

%{
cropSyn = showSyn;
[y x] = find(sum(cropSyn,3)>0);
cropSyn = cropSyn(min(y):max(y),min(x):max(x),:);
image(uint8(cropSyn * 1)),pause(.01)
%}

%{
imwrite(uint8(cropSyn),sprintf('%ssynPos_%s.png',synDir,'907allAxons'))

imwrite(uint8(cropSyn),sprintf('%ssynPos_%06.0f.png',synDir,[ 0 ]))


imwrite(uint8(cropSyn),sprintf('%ssynPosNoSyn_%06.0f.png',synDir,showCell))

%}

