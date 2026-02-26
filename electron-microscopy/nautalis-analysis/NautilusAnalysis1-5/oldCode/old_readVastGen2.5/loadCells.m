
%SPN = 'D:\LGNs1\Segmentation\VAST\S4\joshm\exports\export_fullMip3_14+07+26\'


SPN = GetMyDir;
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
path(path,sprintf('%s\\anaFunc',pwd))

path(path,sprintf('%s',pwd))


%% Read in vast Colors
useSaved = 1;

if ~exist([MPN 'vastSubs.mat'],'file')
  tif2point(SPN,1); %make file for each plane
  stackObs(SPN); %make objects
end

if useSaved & exist([MPN 'obI.mat'],'file')
    load([MPN 'obI.mat'])
else
    obI = makeOBI(SPN,MPN)
end

if ~exist('cellList','var')
    cellList = obI.cell.name;
end


%% Down sample vastSubs

if useSaved & exist([MPN 'dsObj.mat'],'file')
    load([MPN 'dsObj.mat'])
else
   dsDim = [1 1 4];
   dsObj = downSampObj(MPN, dsDim);
end






