
%SPN = 'D:\LGNs1\Segmentation\VAST\S4\joshm\exports\export_fullMip3_14+07+26\'


%SPN = GetMyDir;
SPN = 'D:\otherPeoplesData\Alyssa\VAST_Exports_AlyssaTest1\'
MPN = [SPN(1:end-1) '_mat\']
if ~exist(MPN,'dir'),mkdir(MPN),end

useSaved = 1;


%{

mpSize = matlabpool('size');
if ~mpSize
    'opening matlab pool'
    matlabpool close force
    matlabpool 
end

%}

%% see if folder is stable for 10 sec

oldDir = dir(SPN);

while 1
    pause(10)
    
    newDir = dir(SPN);
    if length(oldDir) == length(newDir)
       'folder stable'
       break
    else
       sprintf('File number increased too %d',length(newDir)) 
    end
    oldDir = newDir;
end



%% Read in vast Colors


if ~exist([MPN 'vastSubs.mat'],'file')
  tif2point(SPN,1); %make file for each plane
  stackObs(SPN); %make objects
end

if useSaved & exist([MPN 'obI.mat'],'file')
    load([MPN 'obI.mat'])
else
    obI = makeOBI(SPN,MPN)
end
% 
% if ~exist('cellList','var')
%     cellList = obI.cell.name;
% end


%% Down sample vastSubs

if useSaved & exist([MPN 'dsObj.mat'],'file')
    load([MPN 'dsObj.mat'])
else
   dsDim = [4 4 4];
   dsObj = downSampObj(MPN, dsDim);
end


%% make golgi library
cellList = obI.cell.name;
%cellList = [10	129	162	170]

colMap = hsv(256);
Dim = 1;
fsize = max(cat(1,dsObj.subs),[],1);


golgiDir = ['D:\LGNs1\Segmentation\VAST\S4\cellLibrary\golgi\'];
golgiDir = [MPN 'golgi1\']
if ~exist(golgiDir,'dir'),mkdir(golgiDir),end

for i = 1:length(cellList)
        showCellNames = cellList(i);

        fileName = sprintf('golgi_%05.0f_dim%d.png',showCellNames,Dim);
 if ~exist([golgiDir fileName],'file')
    col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
    col = col(randperm(size(col,1)),:);
    col = [1 1 1];
    I_partCell = showCellSum(obI,dsObj,showCellNames,col,Dim,fsize);
    I_partCell = 256-I_partCell*3;
    image(uint8(I_partCell))
    disp(showCellNames)
    fileName = sprintf('golgi_%05.0f_dim%d.png',showCellNames,Dim);
    imwrite(uint8(I_partCell),[golgiDir fileName])
    pause(.01)
 end
end


%% make OBJ
cellList = obI.cell.name;
cellList = num2cell(cellList(cellList>=100));
renderOb = 0;
%cellList = [10	129	162	170]
objDir = ['D:\LGNs1\Segmentation\VAST\S4\cellLibrary\objFiles\'];

objDir = [MPN 'obFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end

cellSubs = names2Subs(obI,dsObj,cellList);
downSamp = 4;
cellDat.cellList = cellList;
for i = 1:length(cellSubs)
    sub = cellSubs{i};
    obName = cellList{i};
    if ~ischar(obName),obName = num2str(obName);end
    if ~isempty(sub)
    smallSub = shrinkSub(sub,downSamp);
    tic
    fv = subVolFV(smallSub,[],renderOb);
    fileName = sprintf('%sdSamp%d_%s.obj',objDir,downSamp,obName);
    vertface2obj(fv.vertices,fv.faces,fileName,obName);
    toc
    cellDat(i).subs = sub;
    cellDat(i).fv = fv;
    end
    
end




%% View more cells together
cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

%cellList = ([108  129	109	117	162	131	116	137	130	135	106]);

%



%%faciculation
%{
showCellNames = obI.cell.name;
cellNames = num2cell(intersect(obI.cell.name,[200:300 2000:3000]));



cellNames = cellNames((cellNames>2000) & ( cellNames < 3000));

cellNames = num2cell([201 203 cellNames])
showCellNames = cat(2,num2cell(cellNames), {'boutons'});

showCellNames = num2cell([107 112]);
showCellNames = cellNames
showCellNames ={' '};
showCellNames = {};
%}

fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));

viewProps.viewWindow = double([1 1 1; fsize]);
viewProps.viewWindow = double([1 1 1; fsize(1) fsize(2) fsize(3)]);
viewProps.viewWindow = [minVal; fsize];


%}

%{
%%two pops
pop = obI.cell.name;
pop1 = pop((pop<200)&(pop>100));
pop2 = pop((pop>200) & (pop<300));
showCellNames = num2cell([pop1 pop2]);
col = zeros(length(pop1)+length(pop2),3);
col(1:length(pop1),1) = 1;
col(length(pop1)+1:end,2) = 1;
viewProps.viewWindow = double([1500 1 1; fsize]);

%}

col = colMap(ceil((1:length(dsObj))*256/length(dsObj)),:);

%col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);

viewProps.obI = obI;
viewProps.dsObj = dsObj;
%viewProps.col = col;
viewProps.dim = 3;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;



I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum*1))
imwrite(uint8(I_topSum),[cellPicDir 'smallCol2c.png'])

%%
sumDir = [MPN 'sumImages\relay\'];
%cells2Dir(obI,dsObj,partCells,Dim,sumDir);


%%combind images
I_comb = I_preTargCell;
%imwrite(I_preTargCell,[MPN 'sumImages\I_preTargCell.png'])
%imwrite(I_partCell,[MPN 'sumImages\I_partCell.png'])

I_comb(I_comb == 0) = I_targCellSum(I_comb==0);
I_comb = I_comb + I_partCell/3;
image(I_comb)

I_rgb = cat(3,mean(I_preTargCell,3),mean(I_targCellSum,3),mean(I_partCell,3));
image(uint8(I_rgb))
imwrite(uint8(I_rgb),[MPN 'sumImages\axTarg2\' sprintf('%d.png',targCell)])





%% Skeletonize with shortest paths
    cellList = obI.cell.name;
    
disp(sprintf('Cell List = %s',num2str(cellList)));
TPN = [MPN 'skel\'];
TPNview = [TPN 'view\'];
TPNmat = [TPN 'mat\'];
if ~exist(TPN,'dir'),mkdir(TPN);end
if ~exist(TPNview,'dir'),mkdir(TPNview);end
if ~exist(TPNmat,'dir'),mkdir(TPNmat);end

for i = 1:length(cellList)
        cellTarg = cellList(i);
    disp(sprintf('skeletonizing cell %d (%d of %d)',cellTarg,i,length(cellList)))

    skelFile = sprintf('%s%d.mat',TPNmat,cellTarg);
        viewFile = sprintf('%s%d.png',TPNview,cellTarg);
    if ~exist('skelFile','file')

    
    objectSubs = getCellSubs(obI,dsObj,cellTarg);
    seedSub   = ceil(getSeed(obI,cellTarg)/8);
    
    
    if size(objectSubs,1)>1000
    
    startTime = clock;
    
    pass = 1;
    cellFailed = [];
            cellStruct = subs2arbor(objectSubs,seedSub);

    try
        
    catch err
        err
        pass = 0;
        cellFailed = [cellFailed cellTarg]
    end
    
    if pass
        cellView = cellStruct.sideViews{1};
        skel = cellStruct.skel;
        
        
        save(skelFile,'skel')
        imwrite(cellView,viewFile)
        image(cellView)
        
        cellArbors.cellName(i) = cellTarg;
        cellArbors.arbors(i) = rmfield(cellStruct.arbor,'vox');
    end
    
    
    stopTime = clock;
    end
    end
    
    
    
end

%% Save file

%save([TPN 'cellArbors'],'cellArbors');

%% Create Cell views





