
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
colMap = cat(1,[0 0 0],colMap);

viewProps.dim =3;

%% Set Variables
targCells = [ 108 201  907 903];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% n pops
for uc = 1:6
cladeDir = [MPN  'synPos3\breakdown5\'];
if ~exist(cladeDir,'dir'),mkdir(cladeDir),end;

load('.\data\clade\cladePick_six2.mat')

useCells = [conTo([1:4]).tcrList];


colorTable = [ .2 .2 1; 1 1 0; 1 0 0; 2 2 2; 1 0 1 ; 0 1 0];
%colorTable = [1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1];


col = [];
allMembers = [];
useClade = [uc];
inam = sprintf('tcr%d.png',useClade);
for i = 1:length(useClade)
   members = cladePick.members{useClade(i)};
   members = intersect(useCells,members);
   allMembers = [allMembers members];
   col = cat(1,col,repmat(colorTable(useClade(i),:),[length(members),1]));
    
end

%%Reference lists
showCell = allMembers;

showCellNames = cat(2,num2cell(showCell));


%% get positions
[cellID cellPos] = getList_cellPositions;
usePost = [];
for i = 1:length(allMembers)
    targ = find(cellID == allMembers(i));
    if ~isempty(targ)
        usePost = cat(1,usePost,cellPos(targ,:));
    end
end

groupPos = median(usePost,1);

groupX = round(groupPos(2)/.2);
groupY = round(groupPos(1)/.2);

%% Pick synapses
cellList = obI.cell.name;

prePop{2} = [];
postPop{2} =  [];

prePop{1} = [];
postPop{1} =  [];

prePop{3} = [  ];
postPop{3} =  [];


cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end





fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];
%}


viewProps.viewWindow = double([1100 100 1; fsize(1)-100 fsize(2)-100 fsize(3)]);

medY =  2243 
medX =  1015
viewProps.viewWindow = [medY - 750 medX - 750 1;...
    medY + 750 medX + 750 fsize(3)];

groupX = groupX-viewProps.viewWindow(1,2);
groupY = groupY-viewProps.viewWindow(1,1);



%% Display Variables

viewProps.maxScaleFactor = 3;
viewProps.sumScaleFactor = 3;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;



%% Display Cells

I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum*.5))

if 1

  %%
    I = I_topSum;
    gamma = .5;
    I = I.^gamma * 256/256^gamma;
    SE = strel('disk',2);
    I = imdilate(I,SE);
        I = I * 1;

    
    %I = 255-I;
        I(I<0) = 0;

    keepRat = 1;
    %I = uint8keepRat(I*.5)*keepRat + uint8(I*.5)*(1-keepRat);
    
        I = iBoarder(I,20,1000);
    
        
        rLength = 50;
        rWidth = 10;
        mid = round(size(I,1)/2);
        I(mid-rLength:mid+rLength, mid -rWidth : mid + rWidth,:) = 0;
    I(mid-rWidth:mid+rWidth, mid -rLength : mid + rLength,:) = 0;
    I(mid-rLength:mid+rLength, mid -rWidth : mid + rWidth,1) = 1000;
    I(mid-rWidth:mid+rWidth, mid -rLength : mid + rLength,1) = 1000;
    
    %% Draw center
       I(mid-rLength:mid+rLength, mid -rWidth : mid + rWidth,:) = 10000;
    I(mid-rWidth:mid+rWidth, mid -rLength : mid + rLength,:) = 1000;
    I(mid-rLength:mid+rLength, mid -rWidth : mid + rWidth,1) = 0;
    I(mid-rWidth:mid+rWidth, mid -rLength : mid + rLength,1) = 0;

    %% Draw group median
        I(groupY-rLength:groupY+rLength, groupX -rWidth : groupX + rWidth,:) = 10000;
    I(groupY-rWidth:groupY+rWidth, groupX -rLength : groupX + rLength,:) = 1000;
    I(groupY-rLength:groupY+rLength, groupX -rWidth : groupX + rWidth,1) = 0;
    I(groupY-rWidth:groupY+rWidth, groupX -rLength : groupX + rLength,1) = 0;
    
    
    image(uint8(I))
    %%
    %   imwrite(uint8(I),[cellPicDir 'RGCaxonsVsTCRaxons.png'])
    
     imwrite(uint8(I),[cladeDir inam])
    
end





end











return









%% Draw synapses
%%Ball
ballRad = 4;
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
showSyn(:,:,p) =  synImage;

image(uint8(showSyn*2));


end

maskSyn = repmat(sum(showSyn,3),[1 1 3]);

showSyn(~maskSyn) = I_topSum(~maskSyn)/2;
SE = strel('disk',1);
Idia = imdilate(showSyn,SE);
Idia = imdilate(Idia,SE);
image(uint8(Idia*.1))
% rotSyn = imrotate(showSyn,pi);
% image(uint8(rotSyn)*10)

%% Write images


%{
imwrite(uint8(Idia*.1),[cellPicDir 'sixClades_tcr_dim1.png'])
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

