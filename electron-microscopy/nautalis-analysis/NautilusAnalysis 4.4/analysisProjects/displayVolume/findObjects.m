clear all
%% Load data
load('MPN.mat')
%MPN = 'D:\LGNs1\Export\export_KV_LIN_2018+09+20a\';
%MPN = 'D:\LGNs1\Export\export_KV_LIN_morph_2019+3+7K\';

if ~exist('MPN','var')
    MPN = GetMyDir
end

synDir = [MPN 'synPos3\'];
if ~(exist(synDir,'dir')),mkdir(synDir); end

if ~exist('obI','var') | ~exist('dsObj','var')
    disp('loading')
    load([MPN 'obI.mat']);
    load([MPN 'dsObj.mat']);
end

disp('showing')
dSamp = [8 8 4];
colMap = hsv(256);
colMap = cat(1,[0 0 0],colMap);


showCell = {125};

showCellNames = {125};
viewProps.dim =3;
col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);
if size(col,1) == 1;
    col = [1 1 1];
end

%% Pick synapses
prePop{2} = [];
postPop{2} =  [];

prePop{3} = [];
postPop{3} =  [];

prePop{1} = [];
postPop{1} =  [];


cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'); mkdir(cellPicDir), end

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
mid = [12492, 21454, 5875] %massive input on 30001 225
mid = [13496, 19858, 3743] %glomA

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

[I tops] = showCellsAndMore(viewProps);
image(uint8(I*1)+50)

[ty tx] = find(tops>0);
tz = tops(tops>0);


%% get input
[gx gy] = ginput;
clear allPos
for g = 1:length(gx)
dists = sqrt((ty-gy(g)).^2 + (tx-gx(g)).^2);
targD = find(dists==min(dists),1);



foundDist = dists(targD)
foundY = ty(targD);
foundX = tx(targD);
foundZ = tz(targD);

%xyzPos = [foundX foundY foundZ]
xyzPos = [foundY  foundX foundZ];

allPos(g,:) = xyzPos;

end




%%Translate

allPos2 = allPos;

allPos2(:,3) = fsize(3) - allPos2(:,3);
dim = viewProps.dim;
if dim == 1
    dims = [3 2];
elseif dim == 2
    dims = [3 1];
elseif dim == 3
    dims = [1 2];
end


allPos2 = double(allPos2);
allPos2 = allPos2(:,[1 2 3]) ;
   
 dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
% dSamp = dSamp ./ [4 4 1]; %!!!!!!!!!!!!!!!!!!!!!



allPos2(:,1) = allPos2(:,1)/dSamp(1);
allPos2(:,2) = allPos2(:,2)/dSamp(2);
allPos2(:,3) = allPos2(:,3)/dSamp(3);

allPos2 = round(allPos2);
allPos2(allPos2<1) = 1;


sprintf('%.0f %.0f %.0f',allPos2(end,2),allPos2(end,1),allPos2(end,3))


