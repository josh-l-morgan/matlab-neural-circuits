function[] = makeWOI

%%Make work space of interest
TPN = uigetdir;


global tis

SPN = [glob.datDir 'Analysis\Data\preproc\'];
load([SPN 'ptDat.mat']);
load([SPN 'ROI.mat']);
%load([SPN 'ROI2.mat']);
load([SPN 'SOI.mat']);
load([SPN 'GOI.mat']);
load([SPN 'NOI.mat']);
load([SPN 'MOI.mat']);
load([SPN 'COI.mat']);

isAMC = find(tis.cells.type.typeID == 8);
isVGC = intersect(isAMC,find(tis.cells.type.subTypeID == 1));
vgcCids = tis.cids(isVGC)

%%Load all sms
smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];
clear sms
%roiCid = ptDat(:,3);
runCids = vgcCids;
for i = 1:length(runCids);
    cid = runCids(i);
    %%Get distances between nodes
    disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
    fileName = sprintf('sm_cid%d.mat',cid);
    useSM(i) = 1;
    %sm = load([smDir fileName],'skel2skel','nep');
    load([smDir fileName]);
    sms(i).sm = sm;
end

save([TPN '\WOI.mat'],'SPN','ptDat','ROI','SOI','GOI','NOI',...
    'MOI','COI','vgcCids','sms','-v7.3')
