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
targCells = [109 108 201 907 903 170];

%% Find connectivity
conTo = makeConTo(obI,targCells);

%% pick cells
%%Reference lists
axnet1 =  [1000:1500];
axnet2 = [2000:2031 2038:2300];
axnet3 = [2032:2037];
oneAx = 2040;

showCell = [125];

showCellNames = cat(2,num2cell(showCell));

viewProps.dim =1;
col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:)
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
image(uint8(I*1))

[ty tx] = find(tops>0);
tz = tops(tops>0);


%% get input
[gx gy] = ginput;
clear allPos
for g = 1:length(gx)
dists = sqrt((ty-gy(g)).^2 + (tx-gx(g)).^2);
targD = find(dists==min(dists),1);



foundDist = dists(targD);
foundY = ty(targD);
foundX = tx(targD);
foundZ = tz(targD);

%xyzPos = [foundX foundY foundZ]
xyzPos = [foundX  foundZ foundY].* (obI.em.dsRes * 1000) ./ (obI.em.res .* [4 4 1]);

allPos(g,:) = xyzPos;

end
sprintf('%.0f %.0f %.0f',allPos(1),allPos(2),allPos(3))


