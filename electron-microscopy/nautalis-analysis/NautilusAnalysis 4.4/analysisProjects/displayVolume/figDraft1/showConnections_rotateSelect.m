
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


%% Set Variables
targCells = [ 108 201 903 907 ];


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


showCell = [201 907];%        127         268         266         267];

viewProps.dim =2;

showCellNames = cat(2,num2cell(showCell));

col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);

col = [0 1 0; .7 0 1];

%% two pops
% 
% load('.\data\clade\cladePick_six2.mat')
% useClade = cladePick.members{1};
% allCells = [conTo(:).rgcList];
% 
% 
% pop{1} = intersect(conTo(1).rgcList,useClade);
% pop{2} = intersect(conTo(2).rgcList,useClade);
% pop{3} = setdiff(allCells,useClade);



pop{1} = (conTo(2).rgcList);
pop{2} = (conTo(4).rgcList);
% pop{3} = [ 201 ]
% pop{4} = [  907]

colTable = [0 1 0; .5 .2 2; 0 1 1; 1 0 1];

showCell = unique([pop{:}]);
col = zeros(length(showCell),3);
for p = 1:length(pop)
    getPop = pop{p};
    for i = 1:length(getPop)
        targ = find(showCell == getPop(i));
        col(targ,:) = col(targ,:) + colTable(p,:);
    end
end
showCellNames = num2cell(showCell);


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
% 
% showCellNames = num2cell([ pop1 pop2 pop3]);
% col = zeros(length(pop1)+length(pop2) + length(pop3),3);
% col(1:length(pop1),1) = 1;
% col(length(pop1)+ 1 : length(pop1) + length(pop2),2) = 1;
% col(length(pop1) + length(pop2)+ 1:end,3) = 1;
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

prePop{2} = [1025        2035   127];
postPop{2} =  [1025        2035   127];

prePop{1} = [1025        2035   127];
postPop{1} =  [1025        2035   127];

prePop{3} = [  ];
postPop{3} =  [];


cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end





fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];

viewProps.viewWindow(1,1) = fsize(1)/2;% = [minVal; fsize];

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
viewProps.viewWindow = double([fsize(1)/2 1 1; fsize]);
viewProps.viewWindow = double([1 1 1; fsize]);
viewProps.viewWindow(1,1) = double(fsize(1)/2);

%% Display Variables

viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = 2;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;



%% Display Cells

I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum*.5))


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
    
    %Isize = viewProps.viewWindow(2,dims) - viewProps.viewWindow(1,dims)+1;
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
showSyn(:,:,p) =  synImage(viewProps.viewWindow(1,dims(1)):viewProps.viewWindow(2,dims(1)),...
    viewProps.viewWindow(1,dims(2)):viewProps.viewWindow(2,dims(2)));

image(uint8(showSyn*2));


end

maskSyn = repmat(sum(showSyn,3),[1 1 3]);

showSyn(~maskSyn) = I_topSum(~maskSyn)/2;
image(uint8(showSyn*2))
% rotSyn = imrotate(showSyn,pi);
% image(uint8(rotSyn)*10)


if 0
    %%  Make rotation
    getDeg = [0:10:90];
    clear I_rot;
    movieDir = 'D:\LGNs1\Analysis\movies\mixedClade\rot_axBandDb\';
    mkdir(movieDir)
    saveFile = 1;
    c = 0;
    for d = 1:length(getDeg)
        disp(sprintf('Rotating view %d of %d',d,length(getDeg)));
        viewProps.degRot = getDeg(d);
        I = stereoCellsAndMoreFull(viewProps);
        I = uint8(I.^.7 * 4.5);
        image(I),pause(.01)
        if saveFile
            c = c+1;
            fileName = sprintf('%srot_%05.0f.png',movieDir,c);
            imwrite(I,fileName)
        else
            I_rot{d} = I;
        end
    end
    
    while 1
        for d = 1:length(I_rot)
            image(uint8(I_rot{d})),pause(.3)
        end
        for d = length(I_rot)-1:-1:2
            image(uint8(I_rot{d})),pause(.3)
        end
        for d = 1:length(I_rot)
            image(uint8(I_rot{d})),pause(.3)
        end
        for d = length(I_rot)-1:-1:2
            image(uint8(I_rot{d})),pause(.3)
        end
        
    end
    
    
    
end
%% Write images


%{
imwrite(uint8(showSyn*1),[cellPicDir 'rgcABs_CB.png'])
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

