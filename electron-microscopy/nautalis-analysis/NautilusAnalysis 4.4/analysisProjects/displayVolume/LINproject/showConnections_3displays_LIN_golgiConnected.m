
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
targCells = [108];


%% Find connectivity
conTo = makeConTo(obI,targCells);

%% n popsEach possible pairing between RGCs and TCs in the subnetworks of cell A and B was classified according to whether the RGC and TC had the same seed cell association or a different seed cell association.
allEdges = obI.nameProps.edges
preTarg = preTo(allEdges,125);
postTarg = postTo(allEdges,125);
%%Reference lists



%% All cells
allIDs = unique([allEdges(:,1)' allEdges(:,2)' obI.nameProps.cellNum]);
plot(allIDs)



showList = [preTarg(:,1); postTarg(:,1)];
for i = 1:length(showList)

showCell = [showList(i)];

showCellNames = cat(2,num2cell(showCell));

col = [1 1 1];


fsize = double(max(cat(1,dsObj.subs),[],1))+100;
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];

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

fileName = sprintf('golg_%04.0f.png', showCell);
foldName = 'C:\Users\jlmorgan\Documents\LIN\images\lin125_golgis2\';
if ~exist(foldName,'dir'); mkdir(foldName);end
 imwrite(uint8(I_topSum ),[foldName fileName])

end




























