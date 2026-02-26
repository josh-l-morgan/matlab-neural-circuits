
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
viewProps.dim =1;

%% Set Variables
targCells = [ 907];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% n popsEach possible pairing between RGCs and TCs in the subnetworks of cell A and B was classified according to whether the RGC and TC had the same seed cell association or a different seed cell association. 
load('.\data\clade\cladePick_six2.mat')

useCells = [conTo(:).rgcList];
useRGCs = [conTo(:).rgcList];
useTCRs = [conTo(:).tcrList];


colorTable = [ 10 10 10] 

col = [];
allMembers = [];
useClade = [ 1 3 6];
for i = 1:length(targCells)
   members = targCells(i);
  
   allMembers = [allMembers members];
   col = cat(1,col,repmat(colorTable(useClade(i),:),[length(members),1]));
    
end

%%Reference lists
showCell = [907 useRGCs(1)];

showCellNames = cat(2,num2cell(showCell));

col = [1 0 0; 0 1 1];


%% two pops

%{


pop1 = intersect(conTo(1).tcrList,getList_mixedClade);
pop2 = intersect(conTo(2).tcrList,getList_mixedClade);
pop3 = intersect([conTo([3 4 ]).tcrList],getList_mixedClade);


pop1 = setdiff(conTo(1).rgcList,getList_mixedClade);
pop2 = setdiff(conTo(2).rgcList,getList_mixedClade);
pop3 = setdiff([conTo([3 4 ]).rgcList],getList_mixedClade);

% 
% 
% pop1 = setdiff([conTo([1 2 3 4]).rgcList],getList_mixedClade);
% pop2 = intersect([conTo([1 2 3 4]).rgcList],getList_mixedClade);
% pop3 = setdiff([conTo([1 2 3 4]).rgcList],getList_mixedClade);



% 
% 
% pop1 = setdiff(conTo(1).rgcList,getList_mixedClade);
% pop2 = setdiff(conTo(2).rgcList,getList_mixedClade);
% pop3 = setdiff([conTo([3 4 ]).rgcList],getList_mixedClade);
% 


% 
% pop1 = conTo(1).targ;
% pop3 = conTo(2).targ;
% pop2 = [conTo([3 4 5]).targ];

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
col = zeros(length(pop1)+length(pop2) + length(pop3),3);
col(1:length(pop1),1) = 1;
col(length(pop1)+ 1 : length(pop1) + length(pop2),2) = 1;
col(length(pop1) + length(pop2)+ 1:end,3) = 1;
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

% 
% prePop{7} = useRGCs;
% postPop{7} = useTCRs;
% 
% prePop{7} = cladePick.members{1};
% postPop{7} = [108];
% 
% prePop{8} = cladePick.members{3};
% postPop{8} = [108];

synCol  = [ 1 1 1]*20;
%synCol  = [ .3 .3 1; 0 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 1 0 1; 0 1 0];


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
image(uint8(I_topSum*1))
return

%% Draw synapses
%%Ball
ballRad = 5;
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
    showSyn(:,:,c) =  showSyn(:,:,c) + synImage * synCol(p,c);
end

image(uint8(showSyn*2));


end

maskSyn = repmat(sum(showSyn,3),[1 1 3]);

showSyn(~maskSyn) = I_topSum(~maskSyn)/2;
image(uint8(showSyn*.1))
% rotSyn = imrotate(showSyn,pi);
% image(uint8(rotSyn)*10)

%% Write images


%{
imwrite(uint8(showSyn*.1),[cellPicDir 'syn130.png'])
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

