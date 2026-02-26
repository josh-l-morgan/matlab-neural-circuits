%% Run all Programs up to 1/19/09
%%Sets up sequence of programs to run to find puncta on labeled processes

%%Copyright Josh Morgan and Daniel Kerchensteiner 2005-2009

'Must eliminate repeats in dotfinder'

TPN = GetMyDir; %% Retrieves path for Folder that contains image folder (I)

%% Find out what has been done so far
DoAll.doAll = 0;
DoAll.read = 1;
DoAll.skeleton = 1;
DoAll.dotfinder = 1;  % Set up new status varible
DoAll.ratio = 1;
DoAll.mask = 1;
DoAll.round = 1;
DoAll.smartGuide = 1;
DoAll.group = 1;


if exist([TPN 'Status.mat'])
    load([TPN 'Status.mat']) %Retreive previous progress
else
    Status = DoAll;
end



Status = getVars(Status, 'Which programs would you like to run?');
if Status.doAll
    Status = DoAll;
end

save([TPN 'Status.mat'],'Status')
Settings.Status = Status;

%% Set up necessary Directories
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps

%% Get More Variables

if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
end

if ~isfield(Settings,'ImInfo')
    'Need to collect image info at least once'
    Status.read = 1;
end

if Status.dotfinder

    clear v
    v.blockSize = 200
    v.blockBuffer=30;
    v.thresholdStep = 2;
    v.maxDotSize = 343;
    v.minDotSize=3;
    v.percentBackground = 0.95;
    v.punctaThreshold = 1;
    v.peakCutoff= 0.8;
    v.minFinalDotSize = 3;
    v.roundThreshold = 50;

    v = getVars(v,'Define Dotfinder Variables');
    Settings.dotfinder = v;
    save([TPN 'Settings.mat'], 'Settings')
end

if Status.skeleton

    clear v
    v.minObjSize = 50;
    v.minFillSize = 10;
    v.maxSegLength = 5;

    v=getVars(v , 'Define Skeletonization Variables');
    Settings.skeleton = v;
    save([TPN 'Settings.mat'], 'Settings')
end


clear Settings

%% Check out files
if Status.read
    anaRead(TPN)
    Status.read=0
    save([TPN 'Status.mat'],'Status')
end

%% Run Dot Processing
if Status.skeleton
    'Finding Skel'
    anaSk(TPN)
    Status.skeleton=0;
    save([TPN 'Status.mat'],'Status')
end
if Status.dotfinder
    'Finding Dots'
    dotFinder(TPN)
    %anaDF(TPN)
    Status.dotfinder=0;
    save([TPN 'Status.mat'],'Status')
end
if Status.ratio
    'Ratioing'
    anaRa(TPN)
    Status.ratio=0;
    save([TPN 'Status.mat'],'Status')
end
if Status.mask
    'Masking'
    anaMa(TPN)
    anaCB(TPN)
    %     anaFSc(TPN, DPN) %check for shifts
    Status.mask=0;
    save([TPN 'Status.mat'],'Status')
end

'Running SG once'
anaSG(TPN)

if Status.group
    anaGroup(TPN)
    anaMakeUseOnec(TPN)
else
    anaMakeUseOne(TPN)
end



