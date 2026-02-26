
%SPN = 'D:\LGNs1\Segmentation\VAST\S4\joshm\exports\export_fullMip3_14+07+26\'


%SPN = GetMyDir;
SPN = 'D:\LGNs1\Segmentation\rhoanna\run3combo\'
MPN = [SPN(1:end-1) '_mat\']
if ~exist(MPN,'dir'),mkdir(MPN),end

%{

mpSize = matlabpool('size');
if ~mpSize
    'opening matlab pool'
    matlabpool close force
    matlabpool 
end

%}



%% load paths
path(path,sprintf('%s\\analyzeMorphology\\skeletonizeByShortestDistance',pwd))
path(path,sprintf('%s\\viewCells',pwd))
path(path,sprintf('%s\\analyzeColors',pwd))
path(path,sprintf('%s\\analyzeConnectivity',pwd))
path(path,sprintf('%s\\vast2ob',pwd))
path(path,sprintf('%s\\sub2mesh',pwd))
path(path,sprintf('%s',pwd))



%% Read in vast Colors


if ~exist([MPN 'vastSubs.mat'],'file')
  hd52VastSubs(SPN)
end
% 
% if useSaved & exist([MPN 'obI.mat'],'file')
%     load([MPN 'obI.mat'])
% else
%     obI = makeOBI(SPN,MPN)
% end



%% Down sample vastSubs
if useSaved & exist([MPN 'dsObj.mat'],'file')
    load([MPN 'dsObj.mat'])
else
   dsDim = [4 4 4];
   dsObj = downSampObj(MPN, dsDim);
end

obI.cell.name = 1:length(dsObj);
for i = 1:length(dsObj)
    obI.cell.obIDs{i} = i;
end


%% Pick size Range
clear dsSize
for i = 1:length(dsObj)
   dsSize(i,1) = size(dsObj(i).subs,1);  
end
hist(dsSize)

[sortSizes sorted] = sort(dsSize,'descend')
cellList = sorted;

%% make golgi library

colMap = hsv(256);
Dim = 2;
fsize = max(cat(1,dsObj.subs),[],1);


golgiDir = ['D:\LGNs1\Segmentation\VAST\S4\cellLibrary\golgi\'];
golgiDir = [MPN 'golgi\']
if ~exist(golgiDir,'dir'),mkdir(golgiDir),end

for i = 1:length(cellList)
    showCellNames = cellList(i);
    col = [1 1 1];
    I_partCell = showCellSum(obI,dsObj,showCellNames,col,Dim,fsize);
    I_partCell = 256-I_partCell*3;
    image(uint8(I_partCell))
    disp(showCellNames)
    fileName = sprintf('golgi%08.0f_%05.0f_dim%d.png',sortSizes(i), showCellNames,Dim);
    imwrite(uint8(I_partCell),[golgiDir fileName])
    pause(.01)
end


%% make OBJ
cellList = sorted(7:50);
cellList = num2cell(cellList);
renderOb = 0;
%cellList = [10	129	162	170]
objDir = ['D:\LGNs1\Segmentation\VAST\S4\cellLibrary\objFiles\'];

objDir = [MPN 'obFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end

cellSubs = names2Subs(obI,dsObj,cellList);
downSamp = 1;
%cellDat.cellList = cellList;
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






%% Skeletonize with shortest paths
   
    cellList = sorted(7:50);

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
    %seedSub   = ceil(getSeed(obI,cellTarg)/8);
   
    
    
    startTime = clock;
    
    pass = 1;
    cellFailed = [];
            cellStruct = subs2arbor(objectSubs);

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

%% Save file

%save([TPN 'cellArbors'],'cellArbors');

%% Create Cell views





