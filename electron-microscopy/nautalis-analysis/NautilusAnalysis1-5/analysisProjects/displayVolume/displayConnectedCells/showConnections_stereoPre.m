
clear all
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


colMap = hsv(256);
colMap = cat(1,[0 0 0],colMap);


%% Set Variables
targCells = [ 108 201 109 907 903 170];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% pick cells

%%Reference lists

axnet1 =  [1000:1500];
axnet2 = [2000:2031 2038:2300];
axnet3 = [2032:2037];
oneAx = 2044;
allEdges = obI.nameProps.edges(:,[2 1]);
gotList = getList_giantBoutons(MPN);
linked919 = [903 919 9101 9108 9109 9110]

preAx = 2039
showCell = [1033 1023 1023 1028 1104 1105 1106];
showCell = [108 conTo(1).rgcList];

viewProps.dim =1;

showCellNames = cat(2,num2cell(showCell),{'drb'});

col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);

col(1,:) = 1;
%col = [0 .3 1; 0 1 0];
%col(1:3,:) = .3;%[1 0 0; 0 1 0; 0 0 1; 1 1 1];


%% two pops

%{
popA = conTo(1).tcrList;
popB = conTo(2).tcrList;
popA = popA(popA<1000);
popB = popB(popB<1000);

popA = [conTo(1).rgcList];
popB = [conTo(2).rgcList];

popA = [108];
popB = [201];



pop2 = [intersect(popA,popB)];
pop1 = [setdiff(popA,pop2)];
pop3 = [setdiff(popB,pop2)];

% 
% pop1 = pop1(1:3)
% pop2 = pop2(1:3)
% pop3 = pop3(1:3)

% 
% pop1 = [conTo(1).rgcList];%[9102:9107];
% pop2 = [conTo(2).rgcList 2033 ];%[9101 9108 9109 9110 ];
% pop3 = [conTo(1).rgcList ];%[108 201 907 903];
% 
% % pop1 = [conTo(2).tcrList(17:19)] 
% 
% 
% pop1 = [postTo(allEdges,2008)]';
% pop2 = [postTo(allEdges,2044)]';
% pop3 = [];

showCellNames = num2cell([ pop1 pop2 pop3]);
% 
% col = zeros(length(pop1)+length(pop2) + length(pop3),3);
% col(1:length(pop1),1) = 1;
% col(length(pop1)+ 1 : length(pop1) + length(pop2),2) = 1;
% col(length(pop1) + length(pop2)+ 1:end,3) = 1;

popCol = [1 0 0; 1 1 0; .3 .3 2];


targCol = zeros(length(pop1)+length(pop2) + length(pop3),1);
targCol(1:length(pop1)) = 1;
targCol(length(pop1)+ 1 : length(pop1) + length(pop2)) = 2;
targCol(length(pop1) + length(pop2)+ 1:end) = 3;
col = popCol(targCol,:);


%col = cat(1,ones(length(targCells),3),col); 

%col = colMap(ceil((1:length(dsObj))*256/length(dsObj)),:);

%}



%% labeled CB ?
% pop1 = post2035;
% pop2 = cat(2,conTo(1).tcrList, conTo(2).tcrList);
% cbList = [pop1 pop2];
% cbRad =  [ones(1,length(pop1)) * 8 ones(1,length(pop2)) * 5];
% cbCol =   [ones(1,length(pop1)) * 8 ones(1,length(pop2)) * 5];

%% Pick synapses
cellList = obI.cell.name;

prePop{1} = [];
postPop{1} =  [];

prePop{3} = [];
postPop{3} =  [];

prePop{2} = [];
postPop{2} =  [];


cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end





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

%% Display Variables

viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = 3;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;



%% Display Cells

I_topSum = showCellsAndMore(viewProps);


%viewProps.viewWindow = [1200 800 1; 3000 2600 fsize(3)];
%viewProps.pivot = 

viewProps.degRot = -4;
I_stereo1 = stereoCellsAndMore(viewProps);
image(uint8(I_stereo1*3))

viewProps.degRot = +4;
I_stereo2 = stereoCellsAndMore(viewProps);

I_stereo = cat(2,I_stereo1,I_stereo2);
image(uint8(I_stereo))

I_stereo = cat(2,I_stereo2,I_stereo1);
image(uint8(I_stereo*1.5))

set(gcf,'color','k')

getDeg = [-30:1:30];
parfor d = 1:length(getDeg)
    disp(sprintf('Rotating view %d of %d',d,length(getDeg)));
   viewProps.degRot = getDeg(d);
    I_rot{d} = stereoCellsAndMore(viewProps);
    image(uint8(I_rot{d}*2)),pause(.01)
    
end

%%

while 1
    for d = 1:length(I_rot)
        image(uint8(I_rot{d}*2)),pause(.001)
    end
    
    for d = length(I_rot):-1:1
        image(uint8(I_rot{d}*2)),pause(.001)
    end
end

%% Draw synapses
%%Ball
ballRad = 8;
ball = ones(ballRad*2+1,ballRad*2+1,ballRad*2+1);
[y x z] = ind2sub(size(ball),find(ball));
dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
ball(:) = dists;
ball(ball>ballRad) = 0;
ballSum = sum(ball>0,3);
ballSum = ballSum * 200/max(ballSum(:));

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

image(uint8(showSyn*2))
image(uint8keepRat(showSyn * 2))
image(uint8(showSyn/4) + uint8keepRat(showSyn/10))
image(uint8keepRat(showSyn*.1))

saveI = (uint8(showSyn/4) + uint8keepRat(showSyn/10));

% rotSyn = imrotate(showSyn,pi);
% image(uint8(rotSyn)*10)

%% Write images


%{
imwrite(saveI,[cellPicDir 'pop1pop2Seeds.png'])
imwrite(uint8(I_stereo*.066),[cellPicDir 'pop1pop2RGC_stereo.png'])

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

